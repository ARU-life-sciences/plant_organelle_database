# `oatk` assemblies

These are all `oatk` assemblies with default parameters (-c 90 -k 1001).

## Public data or not?

For our purposes we need public data for EBI to collaborate. 

```bash 
for file in $(find . -maxdepth 1 -type d); do
  NAME=$(echo $file | cut -d/ -f2 | tr '_' ' ');
  echo $NAME;
  /software/team301/edirect/esearch -db assembly -query "${NAME}[Organism]" | /software/team301/edirect/efetch -format docsum | grep "FtpPath_GenBank"
done > check_public.txt

# then check which species are not public (no Ftp data!)
grep -B1 '<FtpPath_GenBank>' check_public.txt | grep -v "<Ftp" | grep -v -f - check_public.txt | grep -v "<Ftp" > not_yet_public_species.txt
grep "^[A-Z]" check_public.txt > all_species.txt
grep -Fxv -f not_yet_public_species.txt all_species.txt > public_species.txt
echo "Chlamydomonas reinhardtii" >> public_species.txt
```

And append out outgroup, *Chlamydomonas reinhardtii*
