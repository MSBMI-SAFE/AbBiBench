
dir="/home/localcolabfold/"
OS="linux"


# navigate to the directory where you want to install localcolabfold
cd $dir

# install linux version
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh

# update to the latest
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/update_${OS}.sh -O update_${OS}.sh
chmod +x update_${OS}.sh
# execute it.
./update_${OS}.sh .

