# Guida per lanciare un’analisi con un manifest diverso

Questa guida descrive come eseguire il workflow di tumor purity con un nuovo gruppo di campioni, utilizzando un manifest diverso da `manifest_95.csv`.

## Ordine delle operazioni

1. Preparare e verificare il database Excel.
2. Generare il nuovo manifest con `workflow/scripts/generate_manifest.py`.
3. Controllare il manifest generato.
4. Modificare il `Snakefile`.
5. Modificare eventualmente `workflow/scripts/purity_comparisonOOP.py`.
6. Eseguire un controllo con Snakemake in modalità `--dry-run`.
7. Lanciare l’analisi completa.
8. Controllare i risultati nella nuova cartella dedicata.

---

## 1. Preparare il database Excel

Il database Excel deve contenere le informazioni necessarie a identificare i campioni e a costruire i percorsi dei file.

Lo script `generate_manifest.py` utilizza attualmente le colonne:

- `ID`
- `RUN`

Il valore di `ID` deve corrispondere al nome del campione utilizzato nei file BAM, VCF, ichorCNA e nei risultati degli altri strumenti.

Esempio:

| ID | RUN |
|---|---|
| sample001 | 123 |
| sample002 | 124 |

Il file Excel utilizzato per il confronto di flow cytometry deve contenere almeno:

- `ID`
- `purity`

Gli ID del file Excel devono coincidere esattamente con i valori della colonna `sample_id` del manifest.

---

## 2. Generare il nuovo manifest

Il manifest viene generato dallo script:

```text
workflow/scripts/generate_manifest.py
```

Lo script:

1. Legge il database Excel.
2. Seleziona il secondo foglio Excel se sono presenti più fogli; altrimenti usa il primo.
3. Legge le colonne `ID` e `RUN`.
4. Costruisce i percorsi dei file BAM, VCF e ichorCNA.
5. Controlla che i file esistano.
6. Inserisce nel manifest solo i campioni per cui tutti i file richiesti sono presenti.
7. Salva il file nella cartella `config/` del progetto.

Il nome generato segue il formato:

```text
manifest_<numero_campioni>.csv
```

Per esempio:

```text
manifest_25.csv
```

### Comando

Eseguire il comando dalla directory principale del progetto:

```bash
python workflow/scripts/generate_manifest.py
```

---

## 3. Controllare il manifest generato

Verificare che il file sia stato creato:

```bash
ls -lh config/manifest_*.csv
head config/manifest_25.csv
```

Il manifest deve contenere almeno queste colonne:

```csv
sample_id,bam,vcf,ichor_res
sample001,/path/sample001.bam,/path/sample001.vcf.gz,/path/sample001.params.txt
sample002,/path/sample002.bam,/path/sample002.vcf.gz,/path/sample002.params.txt
```

Controllare che:

- `sample_id` sia valorizzato;
- ogni `sample_id` sia univoco;
- i file indicati in `bam`, `vcf` e `ichor_res` esistano;
- gli ID non contengano spazi aggiuntivi;
- gli ID coincidano con quelli del file Excel di flow cytometry;
- il numero di campioni sia quello atteso.

Per controllare il numero di righe:

```bash
wc -l config/manifest_25.csv
```

Il numero di campioni è il numero di righe meno la riga di intestazione.

---

## 4. Modificare `generate_manifest.py` se necessario

Se il nuovo database o la struttura delle cartelle è diversa, modificare il file:

```text
workflow/scripts/generate_manifest.py
```

### File Excel di input

Attualmente il file viene cercato in:

```python
excel_file = Path("/mnt/drive02/Tumor_purity_pipeline/db_last.xlsx")
```

Sostituire il percorso se si utilizza un database diverso:

```python
excel_file = Path("/percorso/del/nuovo_database.xlsx")
```

### Nome del foglio Excel

Il codice attuale seleziona automaticamente il secondo foglio quando sono presenti più fogli:

```python
if len(sheet_names) > 1:
    sheet_name = sheet_names[1]
else:
    sheet_name = sheet_names[0]
```

Se necessario, indicare esplicitamente il foglio:

```python
sheet_name = "nome_del_foglio"
df = pd.read_excel(excel_file, sheet_name=sheet_name)
```

### Nomi delle colonne Excel

Se il nuovo database non usa `ID` e `RUN`, modificare:

```python
sample_id = str(row['ID'])
run = str(row['RUN'])
```

Per esempio:

```python
sample_id = str(row['Sample_ID'])
run = str(row['Sequencing_Run'])
```

### Percorsi dei file

Se la struttura delle cartelle è diversa, modificare le variabili:

```python
bam_base_folder = base_dir / 'prunity_share' / 'SEQUENCING' / 'AIO_PANEL' / 'BAM'
vcf_base_folder = base_dir / 'drive02' / 'UMA'
ichor_base_folder = base_dir / 'drive02' / 'UMA'
```

Se necessario, modificare anche i percorsi costruiti per ogni campione:

```python
bam_path = bam_base_folder / f"run_{run}" / f"remDup_recal_{sample_id}.bam"
vcf_path = vcf_base_folder / f"run_{run}" / "Mutect2" / f"mutect2_{sample_id}.vcf.gz"
ichor_path = ichor_base_folder / f"run_{run}" / "ichorCNA" / f"{sample_id}.params.txt"
```

### Cartella di output

Il manifest viene salvato attualmente in:

```python
out_dir = Path("/mnt/drive02/Tumor_purity_pipeline/config/")
output_csv = f"manifest_{len(manifest_df)}.csv"
```

È consigliabile mantenere il file nella cartella `config/` del progetto.

---

## 5. Modificare il `Snakefile`

Supponiamo che il nuovo manifest sia:

```text
config/manifest_25.csv
```

Il file da modificare è:

```text
Snakefile
```

### 5.1 Aggiornare il manifest iniziale

Cambiare:

```python
manifest_df = pd.read_csv("config/manifest_95.csv")
```

in:

```python
manifest_df = pd.read_csv("config/manifest_25.csv")
```

Questa sezione crea la lista `SAMPLES` e le associazioni tra campioni, BAM, VCF e risultati ichorCNA.

### 5.2 Creare una nuova cartella per i risultati

Non utilizzare la cartella dell’analisi precedente, per evitare di sovrascrivere i risultati.

Per esempio, sostituire:

```text
results/purity_comparison_95/
```

con:

```text
results/purity_comparison_25/
```

È possibile usare anche un nome descrittivo, come:

```text
results/purity_comparison_normali/
```

### 5.3 Aggiornare la regola `all`

Esempio:

```python
rule all:
    input:
        "results/purity_comparison_25/purity_comparison_25.csv",
        "results/purity_comparison_25/plots",
        expand("results/ascat_sc/{sample}/ascat.done", sample=SAMPLES),
        expand("results/purecn/{sample}/{sample}.csv", sample=SAMPLES)
```

### 5.4 Aggiornare la regola `purity_comparison`

Cambiare il manifest:

```python
manifest = "config/manifest_25.csv"
```

Aggiornare anche `output` e `visualization_dir`:

```python
output:
    report = "results/purity_comparison_25/purity_comparison_25.csv",
    plots = directory("results/purity_comparison_25/plots")
params:
    flow_db = "db_last.xlsx",
    ascat_dir = "results/ascat_sc/",
    purecn_dir = "results/purecn/",
    visualization_dir = "results/purity_comparison_25/plots"
```

### Nota importante sui risultati intermedi

Se i nuovi campioni hanno gli stessi `sample_id` di una precedente analisi, anche le cartelle:

```text
results/ascat_sc/<sample_id>/
results/purecn/<sample_id>/
```

potrebbero contenere risultati precedenti.

Per mantenere completamente separate due analisi, si possono usare cartelle dedicate, per esempio:

```text
results/ascat_sc_25/
results/purecn_25/
```

In questo caso è necessario aggiornare tutti i percorsi corrispondenti nel `Snakefile` e in `purity_comparisonOOP.py`.

---

## 6. Modificare il file Excel di flow cytometry

Nel `Snakefile`, il file Excel è indicato nella regola `purity_comparison`:

```python
flow_db = "db_last.xlsx"
```

Se si utilizza un file diverso, modificarlo, per esempio:

```python
flow_db = "/mnt/drive02/Tumor_purity_pipeline/flowdb_nuovo_gruppo.xlsx"
```

Il file deve contenere almeno le colonne:

```text
ID
purity
```

Lo script `purity_comparisonOOP.py` legge automaticamente il primo foglio e cerca le colonne `ID` e `purity`, ignorando la differenza tra maiuscole e minuscole.

Se il foglio corretto non è il primo, modificare:

```python
sheet_name = "nome_del_foglio"
df = pd.read_excel(self.flow_db_path, sheet_name=sheet_name)
```

Se le colonne hanno nomi diversi, aggiornare la ricerca delle colonne nello stesso script.

---

## 7. Controllare `purity_comparisonOOP.py`

Il file:

```text
workflow/scripts/purity_comparisonOOP.py
```

non deve normalmente essere modificato se il nuovo manifest e il nuovo Excel rispettano la struttura attesa.

Lo script si aspetta:

### Colonne del manifest

```text
sample_id
bam
vcf
ichor_res
```

La colonna `quality_flag` non è obbligatoria: se non esiste, viene creata automaticamente con valore `PASS`.

### Risultati ASCAT.sc

Per ogni campione viene cercato:

```text
results/ascat_sc/<sample_id>/summary_project_<sample_id>.txt
```

### Risultati PureCN

Per ogni campione viene cercato:

```text
results/purecn/<sample_id>/<sample_id>.csv
```

### Risultati ichorCNA

Il percorso viene letto direttamente dalla colonna `ichor_res` del manifest.

Se i file prodotti da ASCAT.sc o PureCN hanno una struttura diversa, modificare le funzioni:

- `load_ASCAT_results()`
- `load_PureCN_results()`
- `load_ichor_results()`

Controllare in particolare i nomi delle cartelle, dei file e delle colonne contenenti la purezza.

---

## 8. Eseguire un controllo preliminare

Prima dell’analisi completa eseguire un dry-run:

```bash
snakemake --dry-run
```

Oppure:

```bash
snakemake --dry-run --cores 4
```

Il dry-run permette di verificare:

- i campioni individuati;
- i percorsi BAM e VCF;
- le regole che verranno eseguite;
- la nuova cartella dei risultati;
- eventuali file mancanti.

Non procedere con l’analisi completa finché gli errori del dry-run non sono stati risolti.

---

## 9. Lanciare l’analisi completa

Dopo il dry-run, eseguire:

```bash
snakemake --cores 4
```

Il numero di core può essere modificato in base alle risorse disponibili:

```bash
snakemake --cores 8
```

Il workflow eseguirà, nell’ordine:

1. ASCAT.sc per ogni campione.
2. Generazione degli intervalli PureCN.
3. Generazione della coverage dei normali.
4. Generazione della NormalDB.
5. Generazione della coverage tumorale.
6. Esecuzione di PureCN.
7. Confronto tra flow cytometry, ASCAT.sc, ichorCNA e PureCN.
8. Creazione del report finale e dei grafici.

---

## 10. Dove trovare i risultati

Se è stata utilizzata la cartella `results/purity_comparison_25/`, i risultati principali saranno:

```text
results/purity_comparison_25/purity_comparison_25.csv
results/purity_comparison_25/plots/
results/purity_comparison_25/ichorcna_cutoff_summary.csv
results/purity_comparison_25/quality_filter_report.txt
```

### Report finale

Il file principale è:

```text
results/purity_comparison_25/purity_comparison_25.csv
```

Contiene i valori di purezza per campione, tra cui:

- `purity_flow`
- `purity_ascat`
- `purity_ichorCNA`
- `purity_PureCN`

### Grafici

I grafici vengono salvati in:

```text
results/purity_comparison_25/plots/
```

### Report di qualità

Il file:

```text
results/purity_comparison_25/quality_filter_report.txt
```

riporta il numero totale di campioni, i campioni con problemi in PureCN e i campioni presenti nel risultato finale.

### Report del cutoff ichorCNA

Il file:

```text
results/purity_comparison_25/ichorcna_cutoff_summary.csv
```

contiene il riepilogo delle statistiche prima e dopo l’applicazione del cutoff ichorCNA.

### Risultati ASCAT.sc

```text
results/ascat_sc/<sample_id>/
```

Il file principale utilizzato dal confronto è:

```text
summary_project_<sample_id>.txt
```

### Risultati PureCN

```text
results/purecn/<sample_id>/
```

Il file principale utilizzato dal confronto è:

```text
<sample_id>.csv
```

### Log

I log si trovano nella cartella:

```text
logs/
```

In particolare:

```text
logs/ascat_sc/
logs/purecn/
logs/purity_comparison.log
```

---

## 11. Checklist finale

Prima di lanciare il workflow controllare:

- [ ] Il database Excel contiene le colonne corrette.
- [ ] Il foglio Excel selezionato è quello corretto.
- [ ] `generate_manifest.py` utilizza i percorsi corretti.
- [ ] Il nuovo manifest è stato generato nella cartella `config/`.
- [ ] Il manifest contiene `sample_id`, `bam`, `vcf` e `ichor_res`.
- [ ] I file indicati nel manifest esistono.
- [ ] Gli ID del manifest coincidono con quelli del database flow cytometry.
- [ ] Nel `Snakefile` è stato aggiornato il percorso del manifest iniziale.
- [ ] Nella regola `purity_comparison` è stato aggiornato il percorso del manifest.
- [ ] Il file Excel flow cytometry è quello corretto.
- [ ] È stata creata una nuova cartella per i risultati.
- [ ] La nuova cartella è stata aggiornata nella regola `all`.
- [ ] La nuova cartella è stata aggiornata nella sezione `output` di `purity_comparison`.
- [ ] La nuova cartella è stata aggiornata in `visualization_dir`.
- [ ] È stato eseguito `snakemake --dry-run`.
- [ ] I risultati precedenti non verranno sovrascritti.

> **Nota importante:** prima di eseguire il workflow, verificare sempre che i percorsi dei file e i nomi delle colonne siano coerenti tra database Excel, manifest, `Snakefile` e script Python.
