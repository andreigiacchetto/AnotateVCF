import argparse
import csv
import gzip
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Optional, TextIO, Tuple, Union

logger = logging.getLogger(__name__)

# Mapeia (cromossomo_normalizado, posição) -> (rs, snp_name) ou snp_name
SnpEntry = Union[Tuple[str, str], str]
SnpDict = Dict[Tuple[str, str], SnpEntry]

_CHR_PREFIX_RE = re.compile(r"^chr", re.IGNORECASE)


def normalize_chrom(chrom: str) -> str:
    """Normaliza o identificador do cromossomo removendo espaços e prefixo 'chr'."""
    return _CHR_PREFIX_RE.sub("", chrom.strip())


def open_text_file(path: Path, mode: str) -> TextIO:
    """Abre arquivos de texto comuns ou compactados em gzip (.gz) com UTF-8."""
    text_mode = mode if "t" in mode else f"{mode}t"
    if path.suffix.lower() == ".gz":
        return gzip.open(path, text_mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def validate_input_file(path: Path, label: str) -> None:
    """Garante que o arquivo de entrada existe e é um arquivo regular."""
    if not path.is_file():
        raise FileNotFoundError(f"Arquivo {label} não encontrado: {path}")


def process_csv(csv_file_path: Path) -> SnpDict:
    """Lê o arquivo CSV de anotação e retorna um dicionário de SNPs indexado por (cromossomo, posição)."""
    validate_input_file(csv_file_path, "CSV")
    snp_dict: SnpDict = {}
    duplicates = 0
    required_cols = {"chromosome", "position", "rs", "SNP_name"}

    logger.info("Iniciando a leitura do arquivo CSV: %s", csv_file_path)
    with open_text_file(csv_file_path, "r") as csv_file:
        csv_reader = csv.DictReader(csv_file)
        if not csv_reader.fieldnames:
            raise ValueError(f"Arquivo CSV vazio ou sem cabeçalho: {csv_file_path}")

        missing_cols = required_cols - set(csv_reader.fieldnames)
        if missing_cols:
            raise ValueError(
                f"Colunas obrigatórias ausentes no CSV ({', '.join(sorted(missing_cols))}). "
                f"Esperado: {', '.join(sorted(required_cols))}."
            )

        for row in csv_reader:
            chrom = (row.get("chromosome") or "").strip()
            pos = (row.get("position") or "").strip()
            if not chrom or not pos:
                continue

            rs = (row.get("rs") or "").strip()
            snp_name = (row.get("SNP_name") or "").strip()
            key = (normalize_chrom(chrom), pos)

            if key in snp_dict:
                duplicates += 1
            snp_dict[key] = (rs, snp_name)

    if duplicates:
        logger.warning(
            "Foram encontradas %d entradas duplicadas de (cromossomo, posição) no CSV; "
            "a última ocorrência foi mantida.",
            duplicates,
        )
    logger.info("Leitura do CSV concluída. Total de SNPs carregados: %d", len(snp_dict))
    return snp_dict


def process_bpm(bpm_file_path: Path) -> SnpDict:
    """Lê o arquivo BPM (formato TSV com seção [Data]) e retorna um dicionário de SNPs."""
    validate_input_file(bpm_file_path, "BPM")
    snp_dict: SnpDict = {}
    duplicates = 0

    logger.info("Iniciando a leitura do arquivo BPM (TSV): %s", bpm_file_path)
    with open_text_file(bpm_file_path, "r") as bpm_file:
        found_data_section = False
        for line in bpm_file:
            if line.strip().startswith("[Data]"):
                found_data_section = True
                break

        if not found_data_section:
            raise ValueError(f"Seção '[Data]' não encontrada no arquivo BPM: {bpm_file_path}")

        bpm_reader = csv.DictReader(bpm_file, delimiter="\t")
        headers = bpm_reader.fieldnames
        if not headers:
            raise ValueError(f"Cabeçalho ausente após a seção '[Data]' no arquivo BPM: {bpm_file_path}")

        header_map = {key.strip().lower(): key for key in headers}
        chrom_key = header_map.get("chr")
        pos_key = header_map.get("position")
        snp_name_key = header_map.get("snp name")

        if not chrom_key or not pos_key or not snp_name_key:
            raise ValueError(
                "Cabeçalhos 'Chr', 'Position' ou 'SNP Name' não foram encontrados no arquivo BPM."
            )

        for row in bpm_reader:
            first_val = (row.get(headers[0]) or "").strip()
            # Interrompe se iniciar outra seção (ex.: [Controls]) após [Data]
            if first_val.startswith("[") and first_val.endswith("]"):
                break

            chrom = (row.get(chrom_key) or "").strip()
            pos = (row.get(pos_key) or "").strip()
            if not chrom or not pos:
                continue

            snp_name = (row.get(snp_name_key) or "").strip()
            key = (normalize_chrom(chrom), pos)

            if key in snp_dict:
                duplicates += 1
            snp_dict[key] = snp_name

    if duplicates:
        logger.warning(
            "Foram encontradas %d entradas duplicadas de (cromossomo, posição) no BPM; "
            "a última ocorrência foi mantida.",
            duplicates,
        )
    logger.info("Leitura do BPM (TSV) concluída. Total de SNPs carregados: %d", len(snp_dict))
    return snp_dict


def _append_info_field(info: str, tag: str, value: str) -> str:
    """Adiciona ou atualiza uma chave no campo INFO de uma linha VCF."""
    entry = f"{tag}={value}"
    if not info or info == ".":
        return entry
    return f"{info};{entry}"


def process_vcf(
    vcf_file_path: Path,
    output_file_path: Path,
    snp_dict: SnpDict,
    is_csv: bool,
    strict_vcf: bool = False,
) -> None:
    """Anota um arquivo VCF (ou .vcf.gz) com os dados carregados de CSV ou BPM."""
    validate_input_file(vcf_file_path, "VCF de entrada")
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Iniciando a leitura do arquivo VCF: %s", vcf_file_path)
    total_lines = 0
    variant_lines = 0
    substitutions_made = 0
    info_header_written = False

    with open_text_file(vcf_file_path, "r") as vcf_file, open_text_file(
        output_file_path, "w"
    ) as output_file:
        for raw_line in vcf_file:
            total_lines += 1
            line = raw_line.rstrip("\r\n")
            if not line:
                continue

            if line.startswith("#"):
                if line.startswith("##INFO=<ID=SNP_NAME,"):
                    info_header_written = True

                if line.startswith("#CHROM"):
                    if strict_vcf:
                        if not info_header_written:
                            output_file.write(
                                '##INFO=<ID=SNP_NAME,Number=1,Type=String,Description="SNP name from annotation file">\n'
                            )
                        output_file.write(line + "\n")
                    else:
                        header = line.split("\t")
                        header.insert(0, "SNP_name")
                        output_file.write("\t".join(header) + "\n")
                else:
                    output_file.write(line + "\n")
                continue

            variant_lines += 1
            fields = line.split("\t")
            if len(fields) < 3 or (strict_vcf and len(fields) < 8):
                raise ValueError(
                    f"Linha {total_lines} do VCF possui menos colunas do que o esperado ({len(fields)})."
                )

            key = (normalize_chrom(fields[0]), fields[1].strip())
            match: Optional[SnpEntry] = snp_dict.get(key)

            snp_name = ""
            if match is not None:
                if is_csv:
                    rs, snp_name = match  # type: ignore[misc]
                    if rs:
                        fields[2] = rs
                else:
                    snp_name = match  # type: ignore[assignment]
                substitutions_made += 1

            if strict_vcf:
                if snp_name:
                    fields[7] = _append_info_field(fields[7], "SNP_NAME", snp_name)
            else:
                fields.insert(0, snp_name)

            output_file.write("\t".join(fields) + "\n")

            if variant_lines % 10000 == 0:
                logger.info("Processadas %d variantes do arquivo VCF...", variant_lines)

    logger.info("Arquivo VCF modificado salvo como: %s", output_file_path)
    logger.info(
        "Processamento concluído. Linhas totais: %d | Variantes: %d | Anotações aplicadas: %d",
        total_lines,
        variant_lines,
        substitutions_made,
    )


def build_parser() -> argparse.ArgumentParser:
    """Constrói e retorna o parser de argumentos da linha de comando."""
    parser = argparse.ArgumentParser(
        description="Anota um arquivo VCF (ou .vcf.gz) usando um arquivo CSV ou BPM (TSV)."
    )
    source_group = parser.add_mutually_exclusive_group(required=True)
    source_group.add_argument("-csv", help="Caminho para o arquivo CSV", type=Path)
    source_group.add_argument("-bpm", help="Caminho para o arquivo BPM (TSV)", type=Path)

    parser.add_argument(
        "-vcf_input",
        help="Caminho para o arquivo VCF de entrada (.vcf ou .vcf.gz)",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "-vcf_output",
        help="Caminho para o arquivo VCF de saída (.vcf ou .vcf.gz)",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--strict-vcf",
        action="store_true",
        help=(
            "Mantém a especificação VCF padrão adicionando SNP_NAME no campo INFO "
            "em vez de inserir uma nova coluna antes de #CHROM."
        ),
    )

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-q", "--quiet", action="store_true", help="Exibe apenas avisos e erros.")
    verbosity.add_argument("-v", "--verbose", action="store_true", help="Exibe logs detalhados.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    log_level = logging.WARNING if args.quiet else (logging.DEBUG if args.verbose else logging.INFO)
    logging.basicConfig(level=log_level, format="%(levelname)s: %(message)s")

    try:
        if args.csv:
            snp_dict = process_csv(args.csv)
            is_csv = True
        else:
            snp_dict = process_bpm(args.bpm)
            is_csv = False

        process_vcf(
            vcf_file_path=args.vcf_input,
            output_file_path=args.vcf_output,
            snp_dict=snp_dict,
            is_csv=is_csv,
            strict_vcf=args.strict_vcf,
        )
    except (FileNotFoundError, ValueError, OSError) as exc:
        logger.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
