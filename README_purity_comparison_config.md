# Guida rapida per cambiare manifest e database flow cytometry

Se vuoi usare un manifest diverso oppure un file Excel con i dati di flow cytometry diversi, modifica questi punti del workflow.

## 1) Manifest diverso

### File da aggiornare
- [Snakefile](Snakefile)
- [workflow/scripts/purity_comparisonOOP.py](workflow/scripts/purity_comparisonOOP.py)

### Variabili da modificare

Nel file [Snakefile](Snakefile):
- riga 7: cambia il path del manifest usato all'inizio del workflow
  - esempio: `pd.read_csv("config/manifest_normali.csv")`
- riga 244: cambia il path del manifest passato alla regola `purity_comparison`
  - esempio: `manifest = "config/manifest_normali.csv"`

Nel file [workflow/scripts/purity_comparisonOOP.py](workflow/scripts/purity_comparisonOOP.py):
- riga 43: il codice legge il manifest con `pd.read_csv(self.manifest_path)`
- righe 44 e 57-58: controlla che la colonna `sample_id` sia presente e che `quality_flag` esista

### Colonne richieste nel manifest
Il manifest deve contenere almeno queste colonne:
- `sample_id`
- `bam`
- `vcf`
- `ichor_res`

Se il tuo file ha nomi delle colonne diversi, devi aggiornare anche il codice che legge questi valori.

## 2) File Excel di flow cytometry diverso

### File da aggiornare
- [Snakefile](Snakefile)
- [workflow/scripts/purity_comparisonOOP.py](workflow/scripts/purity_comparisonOOP.py)

### Variabili da modificare

Nel file [Snakefile](Snakefile):
- riga 250: cambia il percorso del file Excel di flow cytometry
  - esempio: `flow_db = "/mnt/drive02/Tumor_purity_pipeline/flowdb_normali.xlsx"`

Nel file [workflow/scripts/purity_comparisonOOP.py](workflow/scripts/purity_comparisonOOP.py):
- riga 74: se il file Excel ha un foglio diverso, cambia la lettura con `pd.read_excel(self.flow_db_path, sheet_name=...)`
- riga 77: controlla che gli ID nel file Excel coincidano con i campioni del manifest
- riga 83: il codice seleziona la colonna `purity` e la rinomina in `purity_flow`

### Colonne richieste nel file Excel di flow cytometry
Il file Excel deve contenere almeno queste colonne:
- `ID`
- `purity`

Se il tuo file usa altri nomi per queste colonne, devi aggiornare le righe sopra indicate.

## 3) Esempio minimo di struttura

### Manifest CSV
```csv
sample_id,bam,vcf,ichor_res
sample001,/path/to/sample001.bam,/path/to/sample001.vcf.gz,/path/to/sample001.params.txt
sample002,/path/to/sample002.bam,/path/to/sample002.vcf.gz,/path/to/sample002.params.txt
```

### Excel flow cytometry
```text
ID,purity
sample001,0.65
sample002,0.72
```

## 4) Nota importante
Se cambi nome o posizione del manifest o del file Excel, ricordati anche di aggiornare:
- il path nel comando Snakemake
- il path usato dallo script Python
- i nomi delle colonne attese dal codice
