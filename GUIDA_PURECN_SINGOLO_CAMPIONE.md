# Guida per analizzare un singolo campione con PureCN

Questa guida descrive come eseguire il workflow PureCN del progetto `TumorPurity_UMA` su un solo campione, senza lanciare ASCAT.sc, il confronto finale della purezza o la visualizzazione completa.

> **Nota:** il workflow utilizza comunque un manifest CSV. Per analizzare un solo campione, il manifest deve contenere una sola riga, relativa al campione interessato.

---

## 1. Informazioni necessarie

Per il campione servono:

- un identificativo univoco del campione;
- il file BAM tumorale;
- il file VCF delle varianti somatiche;
- gli intervalli PureCN già configurati oppure i file necessari per generarli;
- il pannello NormalDB PureCN, oppure le risorse per generarlo;
- i percorsi corretti definiti in `config/config.yaml`.

Il manifest deve contenere almeno queste colonne:

```text
sample_id,bam,vcf,ichor_res
```

Per PureCN, la colonna `ichor_res` non viene utilizzata direttamente dal comando PureCN, ma è comunque prevista dal manifest e viene letta dal `Snakefile`.

---

## 2. Creare un manifest con un solo campione

Creare un file, per esempio:

```text
config/manifest_sample001.csv
```

Il contenuto deve essere simile al seguente:

```csv
sample_id,bam,vcf,ichor_res
sample001,/percorso/assoluto/sample001.bam,/percorso/assoluto/sample001.vcf.gz,/percorso/assoluto/sample001.params.txt
```

Per una corretta esecuzione, controllare che:

- `sample_id` corrisponda al nome desiderato per la cartella dei risultati;
- il BAM esista e sia leggibile;
- il VCF esista e sia leggibile;
- i percorsi siano assoluti oppure corretti rispetto alla directory di lavoro;
- il VCF sia compatibile con il genoma e il pannello utilizzati da PureCN;
- il campione sia veramente rappresentato da una sola riga.

Controllare i file con:

```bash
ls -lh /percorso/assoluto/sample001.bam
ls -lh /percorso/assoluto/sample001.vcf.gz
```

Se il BAM è indicizzato, verificare anche il relativo indice:

```bash
ls -lh /percorso/assoluto/sample001.bam.bai
```

---

## 3. Verificare la configurazione

Prima di eseguire il workflow, controllare il file:

```text
config/config.yaml
```

I parametri principali utilizzati dal workflow sono:

```yaml
refGen: /percorso/al/genoma/hs37d5.fa
baits_clean: /percorso/all/probes_targets_hg19_clean.bed
mappability: /percorso/mappabilityAlign100mer.bigWig
normal: /percorso/al/bam/di_riferimento.bam
singularity_image: /percorso/alle/immagini
```

Controllare in particolare:

- `refGen`: riferimento genomico corretto;
- `baits_clean`: file BED degli intervalli del pannello;
- `mappability`: file di mappabilità compatibile con il genoma;
- `normal`: BAM normale utilizzato da ASCAT.sc; non è necessario per PureCN, ma è letto dal `Snakefile`;
- `singularity_image`: cartella contenente l’immagine `purecn_2.10.0-amd64.sif`.

Il percorso dell’immagine viene costruito nel `Snakefile` con:

```python
PURECN_SIF = f"{config['singularity_image']}/purecn_2.10.0-amd64.sif"
```

---

## 4. Modificare temporaneamente il manifest nel `Snakefile`

Nel file `Snakefile`, sostituire il manifest utilizzato nella sezione iniziale:

```python
manifest_df = pd.read_csv("config/manifest_95.csv")
```

con:

```python
manifest_df = pd.read_csv("config/manifest_sample001.csv")
```

Anche nella regola `purity_comparison`, se si utilizza il file originale del workflow, il manifest è indicato separatamente:

```python
manifest = "config/manifest_95.csv"
```

Se non si eseguirà la regola `purity_comparison`, questa seconda modifica non è necessaria per l’analisi PureCN. Tuttavia, è consigliabile aggiornarla comunque per evitare riferimenti a un manifest non pertinente.

---

## 5. Come funziona l’esecuzione PureCN

Il `Snakefile` esegue PureCN in più fasi:

1. Generazione degli intervalli ottimizzati.
2. Calcolo della coverage dei campioni normali.
3. Creazione della NormalDB.
4. Calcolo della coverage tumorale del campione.
5. Esecuzione di PureCN sul campione.

Le regole coinvolte sono:

```text
purecn_intervals
purecn_normalcoverage
purecn_normaldb
purecn_tumor_coverage
purecn_run
```

L’output principale del campione sarà:

```text
results/purecn/sample001/sample001.csv
```

---

## 6. Eseguire solo PureCN per il campione

Dalla directory principale del progetto, dopo aver aggiornato il manifest, eseguire il dry-run indicando direttamente l’output del campione:

```bash
snakemake --dry-run --cores 4 results/purecn/sample001/sample001.csv
```

Se il dry-run è corretto, lanciare:

```bash
snakemake --cores 4 results/purecn/sample001/sample001.csv
```

Il target esplicito `results/purecn/sample001/sample001.csv` permette di eseguire solo le regole necessarie per produrre quel risultato, invece di richiedere tutti i risultati definiti dalla regola `all`.

Se si dispone di più risorse:

```bash
snakemake --cores 8 results/purecn/sample001/sample001.csv
```

---

## 7. Comportamento della regola `purecn_run`

La regola `purecn_run` utilizza:

- la coverage tumorale;
- il VCF del campione;
- gli intervalli PureCN;
- la NormalDB;
- il genoma `hg19`;
- il modello `betabin`;
- la segmentazione `PSCBS`.

In caso di errore o non-convergenza di PureCN, il comando nel `Snakefile` è configurato per creare un file CSV di fallback con valori `NA`:

```csv
Sample.ID,purity,ploidy
sample001,NA,NA
```

Per questo motivo, la presenza del file CSV non garantisce necessariamente che PureCN abbia prodotto una stima valida. È indispensabile controllare il log dell’esecuzione.

---

## 8. Dove trovare i risultati

Per `sample001`, il risultato principale è:

```text
results/purecn/sample001/sample001.csv
```

La coverage tumorale è salvata in:

```text
results/purecn/coverage/sample001/sample001_coverage.txt.gz
```

Il log della coverage tumorale è:

```text
logs/purecn/sample001_coverage.log
```

Il log dell’esecuzione PureCN è:

```text
logs/purecn/sample001_run.log
```

Gli intervalli ottimizzati sono salvati in:

```text
resources/purecn_assets/baits_optimized_hg19.txt
```

La coverage dei normali è salvata in:

```text
resources/purecn_assets/coverage/
```

La NormalDB è salvata in:

```text
resources/purecn_assets/normalDB_panel_aio_hg19.rds
```

---

## 9. Controllare il risultato PureCN

Dopo l’esecuzione, controllare il CSV:

```bash
cat results/purecn/sample001/sample001.csv
```

Il file dovrebbe contenere una colonna `purity` e, normalmente, una colonna relativa alla ploidia.

Per verificare la colonna di purezza con Python:

```bash
python - <<'PY'
import pandas as pd

path = "results/purecn/sample001/sample001.csv"
df = pd.read_csv(path)
print(df.columns.tolist())
print(df[['purity']].head())
PY
```

Se il valore di `purity` è `NA`, controllare il log:

```bash
less logs/purecn/sample001_run.log
```

Cercare messaggi relativi a:

- non-convergenza;
- assenza di varianti somatiche sufficienti;
- problemi con PSCBS;
- file VCF non leggibile;
- problemi con la NormalDB;
- intervalli incompatibili;
- errori di memoria o di risorse.

---

## 10. Creare una cartella separata per l’analisi singola

Se il campione è già stato analizzato in precedenza, il workflow potrebbe considerare alcuni file già aggiornati e non rieseguire tutte le regole.

Per evitare di mescolare risultati di analisi diverse, è consigliabile usare una directory di lavoro o di risultati dedicata, per esempio:

```text
results_single_sample/purecn/sample001/
```

Per utilizzare una cartella diversa, occorre modificare i percorsi di output nel `Snakefile`, almeno nelle regole:

- `purecn_tumor_coverage`;
- `purecn_run`.

È necessario modificare anche il target passato a Snakemake. Per esempio:

```bash
snakemake --cores 4 results_single_sample/purecn/sample001/sample001.csv
```

> **Attenzione:** non è sufficiente modificare soltanto il comando finale. I percorsi definiti negli output delle regole devono essere coerenti con il target richiesto.

---

## 11. Esecuzione senza modificare permanentemente il `Snakefile`

Il `Snakefile` attuale definisce il manifest direttamente nel codice. Per eseguire più campioni singolarmente senza dimenticare di ripristinare il file, una possibilità semplice è:

1. creare un manifest per il campione;
2. modificare temporaneamente il percorso del manifest nel `Snakefile`;
3. eseguire il comando Snakemake con il target del campione;
4. ripristinare il percorso originale dopo l’analisi.

Prima di modificare il file, creare una copia:

```bash
cp Snakefile Snakefile.backup
```

Dopo l’esecuzione, ripristinarlo se necessario:

```bash
cp Snakefile.backup Snakefile
```

---

## 12. Checklist finale

Prima di lanciare PureCN su un singolo campione:

- [ ] Il manifest contiene una sola riga di campione.
- [ ] Il manifest contiene `sample_id`, `bam`, `vcf` e `ichor_res`.
- [ ] Il BAM esiste ed è leggibile.
- [ ] Il BAM ha il relativo indice, se richiesto dall’ambiente.
- [ ] Il VCF esiste ed è leggibile.
- [ ] Il genoma di riferimento è quello corretto.
- [ ] Il BED del pannello è corretto.
- [ ] Il file di mappabilità è compatibile.
- [ ] La NormalDB è disponibile oppure può essere generata.
- [ ] L’immagine Singularity PureCN esiste.
- [ ] Il manifest nel `Snakefile` è quello del campione.
- [ ] È stato eseguito il dry-run.
- [ ] È stato usato il target specifico `results/purecn/<sample_id>/<sample_id>.csv`.
- [ ] È stato controllato il log `logs/purecn/<sample_id>_run.log`.
- [ ] Il valore di `purity` non è `NA`.
- [ ] I risultati sono stati salvati in una cartella separata se esiste una precedente analisi dello stesso campione.

---

## Comando minimo riassuntivo

Per un campione chiamato `sample001`:

```bash
# 1. Controllo del workflow
snakemake --dry-run --cores 4 results/purecn/sample001/sample001.csv

# 2. Esecuzione
snakemake --cores 4 results/purecn/sample001/sample001.csv

# 3. Controllo del risultato
cat results/purecn/sample001/sample001.csv

# 4. Controllo del log
less logs/purecn/sample001_run.log
```
