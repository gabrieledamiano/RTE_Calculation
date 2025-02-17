# Modello di Trasferimento Radiativo Atmosferico

Questo repository contiene implementazioni in MATLAB e Python di un modello di trasferimento radiativo atmosferico che calcola la radianza atmosferica utilizzando dati spettroscopici HITRAN. Il tutto finalizzato alla realizzazione del progetto per il corso di Metodi e Tecniche per l'Osservazione della Terra (LM Ingegneria Informatica e delle tecnologie dell'informazione UNIBAS ), a cura del dott. Giuliano Liuzzi.

## Panoramica

Il progetto include due implementazioni principali:
- Uno script MATLAB per l'esecuzione da riga di comando
- Un'applicazione GUI Python con capacità di visualizzazione interattiva

Entrambe le implementazioni eseguono gli stessi calcoli di base:
- Caricamento dei dati atmosferici
- Calcolo delle altezze degli strati
- Calcolo delle densità colonnari dei gas
- Esecuzione dei calcoli di trasferimento radiativo
- Visualizzazione dei risultati inclusi radianza ascendente/discendente e contributi solari

## Prerequisiti

### Per l'implementazione MATLAB
- MATLAB (testato su R2020b o versioni successive)
- File richiesti:
  - `input_atmosfera.mat` e `atmos_data4.mat` (dati atmosferici)
  - File dati HITRAN:
    - `H2O_hitran.txt`
    - `CO2_hitran.txt`
    - `O3_hitran.txt`
    - `NH3_hitran.txt`
    - `HNO3_hitran.txt`

### Per l'implementazione Python
- Python 3.8 o versioni successive
- Pacchetti richiesti:
  ```
  numpy>=1.24.3
  matplotlib>=3.7.1
  scipy>=1.10.0
  tkinter (generalmente incluso in Python)
  ```
- File richiesti (stessi della versione MATLAB):
  - `input_atmosfera.mat` e `atmos_data4.mat`
  - Tutti i file dati HITRAN

## Installazione

1. Installare le dipendenze Python (se si usa la versione Python):
   ```bash
   pip install numpy matplotlib scipy
   ```

2. Posizionare i file dati richiesti nella stessa directory degli script:
   - File dati atmosferici (`input_atmosfera.mat` e `atmos_data4.mat`)
   - File dati HITRAN (tutti i file .txt)

## Utilizzo

### Versione MATLAB
1. Aprire MATLAB
2. Navigare alla directory del progetto
3. Eseguire lo script:
   ```matlab
   >> run_atmosphere_model
   ```

### Versione GUI Python
1. Aprire un terminale/prompt dei comandi
2. Navigare alla directory del progetto
3. Eseguire lo script Python:
   ```bash
   python atmospheric_gui.py
   ```
4. Utilizzare la GUI per:
   - Caricare il file dati atmosferici
   - Calcolare la radianza
   - Visualizzare grafici interattivi
   - Alternare tra temi chiaro/scuro

## Caratteristiche

- Calcolo delle altezze degli strati atmosferici
- Calcolo della densità colonnare dei gas
- Modellazione del trasferimento radiativo includendo:
  - Radianza ascendente
  - Radianza discendente
  - Contributo solare
- Visualizzazione interattiva (versione Python)
- Supporto temi chiaro/scuro (versione Python)


## Funzioni Incluse

Entrambe le versioni includono implementazioni di:
- `blackbodyn`: Calcolo della radiazione del corpo nero
- `lorentzian`: Calcolo del profilo lorentziano
- `calcola_altezze_strati`: Calcolo delle altezze degli strati atmosferici
- `calcola_medie`: Calcolo delle medie degli strati
- `calcola_gas_col`: Calcolo della densità colonnare dei gas

## Output

I programmi generano diverse visualizzazioni:
1. Componenti radiative dirette
2. Componenti riflesse
3. Radianza totale
4. Spettri di assorbimento dei gas
5. Profili di Pressione vs Temperatura
6. Profili dei Mixing Ratio dei gas

## Contributi

Gabriele Damiano e 
Donato Gabriele Carretta

