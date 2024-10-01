# change to a directory where you will be working from
# base albertine or .. above albertine recommended

# download admixture
wget https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz

# unzip folder
tar -xvf admixture_linux-1.3.0.tar.gz

# remove file
rm admixture_linux-1.3.0.tar.gz

cd dist/admixture_linux-1.3.0/

# test to see if commands come up
./admixture

# can reference as shortcut
# MUST DO EVERY SESSION / SHELL SCRIPT
admixture=/work/cooperlab/cyanolyca/dist/admixture_linux-1.3.0/admixture

# should call up program
$admixture

# need to look more into permanently referencing program