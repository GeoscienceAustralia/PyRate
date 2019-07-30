echo "This script will install packages required by PyRate. Continue?" 
    select yn in "Yes" "No"; do
    case $yn in
        Yes ) break;;
        No ) exit;; 
    esac
done

# OS package requirements for Ubuntu 18.04
sudo apt-get update
sudo apt-get -y install \
    gdal-bin \
    libgdal-dev \
    openmpi-bin \
    libopenmpi-dev 

