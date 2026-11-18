# Anotação de Arquivos VCF (`anotate_vcf.py`)

Script em Python para anotar arquivos VCF (inclusive compactados em `.vcf.gz`) com nomes de SNPs (`SNP_name`) e, opcionalmente, identificadores `rs`, a partir de um arquivo CSV ou de um arquivo BPM (formato TSV).

## Requisitos

- Python 3.x
- Utiliza apenas bibliotecas padrão do Python (`argparse`, `csv`, `gzip`, `logging`, `pathlib`, `re`, `sys`, `typing`), sem necessidade de instalar dependências externas.

## Argumentos da Linha de Comando

| Argumento | Obrigatório | Descrição |
| :--- | :--- | :--- |
| `-vcf_input` | Sim | Caminho para o arquivo VCF de entrada (`.vcf` ou `.vcf.gz`). |
| `-vcf_output` | Sim | Caminho para o arquivo VCF de saída (`.vcf` ou `.vcf.gz`). |
| `-csv` | Condicional* | Caminho para o arquivo CSV de anotação. |
| `-bpm` | Condicional* | Caminho para o arquivo BPM (TSV) de anotação. |
| `--strict-vcf` | Não | Mantém a especificação VCF padrão gravando `SNP_NAME` no campo `INFO` em vez de inserir uma nova coluna antes de `#CHROM`. |
| `-q`, `--quiet` | Não | Exibe apenas avisos e erros no log. |
| `-v`, `--verbose` | Não | Exibe logs detalhados. |

> \* `-csv` e `-bpm` são mutuamente exclusivos e exatamente um deles deve ser informado.

## Formatos de Entrada Esperados

### 1. Arquivo CSV (`-csv`)
Arquivo separado por vírgulas (UTF-8) contendo obrigatoriamente as seguintes colunas:
- `chromosome`
- `position`
- `rs`
- `SNP_name`

### 2. Arquivo BPM em TSV (`-bpm`)
Arquivo separado por tabulações (`\t`, UTF-8) contendo uma linha de seção `[Data]` seguida pelo cabeçalho com as colunas (a verificação ignora maiúsculas/minúsculas):
- `Chr`
- `Position`
- `SNP Name`

Linhas anteriores à seção `[Data]` e eventuais seções posteriores (ex.: `[Controls]`) são ignoradas.

### 3. Arquivo VCF (`-vcf_input`)
Arquivo `.vcf` ou `.vcf.gz`. A correspondência é feita por `(cromossomo, posição)` com normalização automática do prefixo `chr` (ou seja, `chr1` e `1` são tratados como equivalentes).

## Exemplos de Uso

Anotando um VCF usando CSV (modo padrão com coluna `SNP_name` no início):

```bash
python anotate_vcf.py -csv snps.csv -vcf_input entrada.vcf -vcf_output saida_anotada.vcf
```

Anotando um `.vcf.gz` usando BPM e preservando o formato VCF padrão (`--strict-vcf`):

```bash
python anotate_vcf.py -bpm manifest.bpm.tsv -vcf_input entrada.vcf.gz -vcf_output saida_anotada.vcf.gz --strict-vcf
```

## Comportamento da Saída

### Modo Padrão
- As linhas de metadados (`##`) são mantidas inalteradas.
- Na linha `#CHROM` e em todas as linhas de variantes, o script insere a coluna `SNP_name` **na primeira posição** (antes de `#CHROM`), deixando-a vazia (`""`) quando não há correspondência.
- No modo `-csv`, quando há correspondência e `rs` não está vazio, a coluna original `ID` também é substituída pelo valor de `rs`.

> **Nota:** Como o modo padrão insere uma coluna antes de `#CHROM`, a saída tem formato tabular personalizado e pode não ser aceita por ferramentas que exigem estritamente a especificação VCF (como `bcftools` ou `vcftools`). Para manter compatibilidade total com essas ferramentas, utilize `--strict-vcf`.

### Modo `--strict-vcf`
- Adiciona a definição `##INFO=<ID=SNP_NAME,Number=1,Type=String,Description="SNP name from annotation file">` ao cabeçalho (caso ainda não exista).
- Mantém todas as colunas originais do VCF intactas e acrescenta `SNP_NAME=<valor>` ao campo `INFO` das variantes encontradas.
- No modo `-csv`, também atualiza a coluna `ID` com `rs` quando disponível.
