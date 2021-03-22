import glob
import os

def main():
    all_files = sorted(glob.glob('/Users/fatemah/Documents/LOPEZ_RESEARCH-pycharm/active_search/iteration*nm*.txt'))
    print(all_files)
    smiles_and_nm = {}
    for file in all_files:
        with open(file, "r") as b:
            print(file)
            for line in b:
                mol = line.split(",")
                smiles_and_nm[mol[0]] = mol[2]
    with open("collected-smiles-data.csv", "w+") as f:
        f.write('smiles,nm\n')
        for key,value in smiles_and_nm.items():
            f.write(key + "," + value)

main()