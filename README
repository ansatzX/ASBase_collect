# ASBase collect 

This repo store scripts which can grab web database to h5 databse.

## usage

1. download data from web [ASBase](http://119.91.135.188:8080)

```
python get_pages.py
```

2. parser html resource and store them to h5py database


```
python parser_asbase.py
```

3. extract infomation from h5 database and construct it to 3D molecule files

```
mkdir pdbs
mkdir xyzs
python extract_3D.py
```

4. transfer 3D file to gaussian input 

first use openbabel convert .pdb to .xyz
```
obabel -ixyz -opdb input.pdb  -O input.xyz
```

then call [geninput](https://github.com/ansatzX/graduation-thesis-BSC/blob/main/tocomp/mole/geninput), a shell script which convert .xyz to gaussian input .gjf

```
geninput input.xyz
```

5. batch computation by PBS or SLURM

```
qsub xxxx
```