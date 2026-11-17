# Odpornost petvhodnih majoritetnih vrat QCA na izpad in premik celic

Seminarska naloga pri predmetu Kvantni celični avtomati in kvantno procesiranje (2026/2027).

Analiziramo napakam odporna petvhodna majoritetna vrata PftMV5 iz članka
O. Selvi in sod., Computers 2026, 15(5), 316, doi: [10.3390/computers15050316](https://doi.org/10.3390/computers15050316).
Vsako okvarjeno različico simuliramo z bistabilnim modelom QCADesigner-E (izčrpno preverjanje vseh 32 vhodnih
kombinacij) in jo primerjamo z brezhibnimi vrati:

- **izpad celice:** odstranimo vsako od 17 notranjih celic (17 simulacij);
- **premik celice:** vsako od 23 celic premaknemo za 1, 2 in 5 nm gor, dol, levo in desno (276 simulacij).

## Datoteke

| Pot | Vsebina |
|---|---|
| `MV5` | izvirna zasnova vrat (QCADesigner) |
| `cells.csv` | oznake celic iz članka ↔ vrstni red celic v datoteki, vloga, položaj, urina cona |
| `engine/main_qca_sim.c` | ukazni simulator (izhaja iz `main_batch_sim.c` iz QCADesigner-E) |
| `engine/build.sh` | prenese izvorno kodo QCADesigner-E in prevede `engine/bin/qca_sim` |
| `sweep.py` | generira okvarjene različice, jih simulira in zapiše `manifest.json` |
| `analyze.py` | iz rezultatov izračuna `segments.csv`, `summary.csv` in `report.md` |
| `reproduce.ps1` | ponovi vse poskuse (semena 1, 2 in 3) |
| `plot.py` | izriše grafa za poročilo (PDF) iz `summary.csv` |

Izvorna koda QCADesigner-E ni v repozitoriju. `engine/build.sh` jo ob prvem zagonu prenese z
[GitHuba](https://github.com/FSillT/QCADesigner-E) pri commitu `3bc90302` v `engine/QCADesigner-E/`
in jo prevede skupaj z `engine/main_qca_sim.c`. Njenih datotek ne spreminjamo.

## Zahteve

- Windows z WSL 2 in distribucijo `Ubuntu` (preizkušeno na Ubuntu 24.04, gcc 13.3). V WSL:
  ```bash
  sudo apt install git build-essential pkg-config libglib2.0-dev
  ```
- Python 3 v sistemu Windows; simulacije in analiza potrebujejo samo standardno knjižnico, izris slik pa `matplotlib`.

## Ponovitev rezultatov

V PowerShellu v korenu repozitorija:

```powershell
powershell -ExecutionPolicy Bypass -File reproduce.ps1
```

Skripta prevede simulator in izvede vse poskuse s semeni 1, 2 in 3 (približno 4 minute na seme).
Posamezno:

```powershell
python sweep.py all --seed 1 --rebuild    # removal | shift | all
python analyze.py runs\<mapa_zagona>      # ponovna analiza zagona
```

Grafa za poročilo (postavitev vrat in valovne oblike so posnetki zaslona iz QCADesigner-E):

```powershell
python plot.py --summary runs\<zagon>\summary.csv --out img
```

Vsak zagon se zapiše v novo mapo `runs\<datum>_<čas>_<vrsta>[_<ime>]\`:
`report.md` (poročilo), `summary.csv` (ena vrstica na različico), `segments.csv` (ena vrstica na različico in
vhodno kombinacijo), `manifest.json` (parametri in zgoščene vrednosti vhodnih datotek), `raw\` (izhod simulatorja)
in `designs\` (okvarjene različice, ki jih lahko odpremo v QCADesigner-E). Rezultati so deterministični.

## Licenca

QCADesigner-E je pod licenco GNU GPL; `engine/main_qca_sim.c` izhaja iz te kode in je pod isto licenco.
