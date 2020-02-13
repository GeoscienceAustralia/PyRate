#   This script is part of the PyRate software package.
#
#   Copyright 2020 Geoscience Australia
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
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

