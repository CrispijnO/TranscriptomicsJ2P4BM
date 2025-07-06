# Genexpressie-analyse duidt op duidelijke verschillen in patiënten met reumatoïde artritis.
<p align=center>
  <img src="Assets/ReumaTitelPlaatje.jpg" width=700px>
  <br>
  Bron: <a href="https://livit.nl/hulp-bij/hand/reuma-in-de-handen/">Livit</a>
</p>
<br>
Door: Crispijn Oppenhuizen<br>
Klas: LBM2-B

# Inhoud
- `Assets/` bestanden voor de opmaak van README.md
- `Bronnen/` gebruikte bronnen
- `Data stewardship/` beschrijving over data beheren op Github
- `Data/` gebruikte data
- `Data/Data_RA_raw` ruwe dataset die gebruikt is voor de analyse
- `Data/aligned Data` aligned data die is geproduceerd van de ruwe dataset
- `Resultaten/` geproduceerde grafieken en tabellen door middel van script
- `Scripts/` gebruikte script voor de analyse

---

# Inleiding
Reumatoïde artritis(RA) is een chronische ontstekingsziekte in de gewrichten, het is een auto-immuun ziekte waarvan de oorzaak nog niet duidelijk is. De ziekte kan uiteindelijk zorgen voor schade aan de gewrichten en organen ([Radu & Bungau, 2021](Bronnen/RaduBungau2021.pdf)). RA is te onderscheiden van andere auto-immuun ziektes door synoviale ontsteking, productie van antilichamen zoals reumafactor(RF) en anti-citrullinated protein antibodies(ACPA) en botvervorming ([Jang et al., 2022](Bronnen/Jang2022.pdf)). De ziekte komt vooral voor in vrouwen en oudere mensen, maar niet exclusief tot deze groepen. De oorzaak van de ziekte is nog niet bekend, maar wel is duidelijk omgevingsfactoren, zoals roken een rol spelen in de ontwikkeling van RA. Ook spelen genetische factoren een grote rol in RA ([Scott et al., 2010](Bronnen/Scott2010.pdf)). Een gen familie in het bijzonder komt voortdurend voor bij patiënten met RA, dit zijn de HLA DRB1 allelen in het MHC op chromosoom 6. De HLA regio codeert voor de MHC klasse 1 en 2 receptoren. Waaronder de alpha en beta ketens van de MHC klasse 2 receptoren. De receptoren zitten op antigeen presenterende cellen en presenteren deze aan T-cellen ([Dedmon, 2020](Dedmon2020.pdf)). Verder is het niet helemaal duidelijk welke genetische factoren er meer een rol spelen in RA. Dit onderzoek probeert hier meer duidelijkheid over te krijgen door de synoviale vloeistof van verschillende mensen te analyseren en hier de verschillen in genexpressie duidelijk te maken.  

# Methode
## Verkregen data
De data dat is gebruikt in dit onderzoek is verkregen door Ilumina sequencing. De sequencing is gedaan op het synoviale vloeistof uit 8 vrouwelijke personen, waarvan 4 een diagnose van RA hebben voor langer dan 12 maanden en 4 zonder RA, waarvan de personen met RA positief testen op ACPA en de gezonden negatief. Verder staan er in tabel 1 meer data over de personen waarvan de samples zijn. 


*Tabel 1 Meta data van personen voor onderzoek*
| Monster ID | Leeftijd | Geslacht   | Status         |
| ---------- | :------: | ---------- | :------------: |
| SRR4785819 | 31       | Vrouwelijk | Controle       |
| SRR4785820 | 15       | Vrouwelijk | Controle       |
| SRR4785828 | 31       | Vrouwelijk | Controle       |
| SRR4785831 | 42       | Vrouwelijk | Controle       |
| SRR4785979 | 54       | Vrouwelijk | Reuma Artritis |
| SRR4785980 | 66       | Vrouwelijk | Reuma Artritis |
| SRR4785986 | 60       | Vrouwelijk | Reuma Artritis |
| SRR4785988 | 59       | Vrouwelijk | Reuma Artritis |

##
