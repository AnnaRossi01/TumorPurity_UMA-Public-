# Copy Number Burden Workflow

## Obiettivo

Questa analisi valuta l'instabilità genomica di ciascun campione tumorale attraverso il Copy Number Burden (CNB), con l'obiettivo di comprendere e giustificare il comportamento degli algoritmi di stima della purezza tumorale su una coorte di pazienti affetti da Mieloma Multiplo.

## Workflow definitivo

### Fase 1: Arbitro esterno (ortogonalità)

I file .cns ottenuti da CNVkit vengono convertiti in formato .seg. Questo passaggio consente di ottenere una rappresentazione indipendente dello stato genomico reale del campione, senza introdurre bias legati agli strumenti oggetto di confronto.

### Fase 2: Classificazione Flat vs Altered

Lo script Python elaborare i file .seg e calcola il Copy Number Burden, definito come la percentuale di genoma alterata all'interno del pannello target. Se il burden è inferiore a una soglia prefissata, il campione viene classificato come "Flat"; se è superiore, come "Altered".

### Fase 3: Stratificazione

Nel database principale contenente i campioni tumorali, viene aggiunta una colonna con la classificazione "Flat" o "Altered" per ciascun ID.

### Fase 4: Analisi critica

Il database viene diviso nei due sottogruppi e le visualizzazioni principali vengono rieseguite separatamente per la coorte Flat e per quella Altered, includendo matrici di correlazione, MAE e RMSE.

---

## Appunti per la tesi: Analisi del Copy Number Burden (CNB)

### 1. Rationale

L'analisi del CNB è stata introdotta per indagare il collasso delle performance degli algoritmi di stima della purezza tumorale (ASCAT.sc, IchorCNA, PureCN) nella coorte di pazienti con Mieloma Multiplo.

Gli strumenti basati sul Copy Number non quantificano direttamente le cellule tumorali, ma stimano la loro frazione analizzando l'instabilità genomica. L'ipotesi centrale è che il loro fallimento non dipenda esclusivamente da errori matematici, ma anche da limiti biologici e tecnici.

- Limite biologico: il MM è eterogeneo e non tutti i pazienti mostrano un genoma fortemente instabile.
- Limite tecnico: il sequenziamento targeted copre una porzione molto ridotta del genoma (0.46 MB). Se le alterazioni non ricadono nelle regioni targettizzate, il pannello appare virtualmente piatto agli algoritmi.

### 2. Metodologia e strumento ortogonale

Per evitare un bias di validazione interna, non sono stati utilizzati i file di segmentazione prodotti dagli strumenti sotto esame. È stato scelto CNVkit come approccio ortogonale e indipendente per generare i file .cns, successivamente convertiti nel formato .seg del Broad Institute.

### 3. Parametri valutati

Lo script Python ha estratto dai file .seg di CNVkit i seguenti parametri:

- seg.mean (Log2 Ratio): misura la deviazione del numero di copie rispetto al baseline diploide.
- num.mark: indica il numero di marcatori o sonde inclusi in ciascun segmento.
- Copy Number Burden (CNB): frazione o percentuale del pannello target interessata da alterazioni del numero di copie.

Formula:

$$
CNB = \frac{\text{somma dei num.mark dei segmenti alterati}}{\text{somma totale dei num.mark del pannello}}
$$

### 4. Cut-off e soglie utilizzate

Sono state definite due soglie analitiche:

- Soglia di alterazione del segmento: valori di Log2 Ratio superiori a 0.2 o inferiori a -0.2.
- Soglia di classificazione del genoma: CNB pari al 15% (0.15).

Classificazione:

- Flat: CNB < 0.15
- Altered: CNB >= 0.15

### 5. Obiettivo della stratificazione

La divisione in due gruppi consente di ricalcolare metriche di validazione come Spearman correlation, MAE e RMSE separatamente per i campioni Flat e Altered.

- Se le performance migliorano nel gruppo Altered, l'utilità degli algoritmi risulta confermata, ma limitata ai pazienti con alto burden genomico.
- Se le performance rimangono inaccettabili in entrambi i gruppi, ciò suggerisce un limite strutturale insormontabile del sequenziamento targeted a 0.46 MB.
