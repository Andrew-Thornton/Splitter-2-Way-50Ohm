# Splitter-2-Way-50Ohm

This GitHub repository provides a design for a **50-ohm RF splitter**.  
It currently uses **KiCad 10.0** and is being tested on **Ubuntu 24.04**.


---


## Installation Instructions

### 1. General Git Repository Setup

Initialize the repository and its submodules:

```bash
git submodule update --init --recursive
```

### 2. Kicad
Add the KiCad PPA, update your package list, and install KiCad:
```bash
sudo add-apt-repository ppa:kicad/kicad-10.0-releases  
sudo apt update  
sudo apt install kicad  
```

### 3. OpenEMS
Install Required Packages:
```bash
sudo apt-get update  
sudo apt-get install build-essential \  
                       cmake \  
                       git \  
                       libhdf5-dev \  
                       libboost-all-dev \  
                       libcgal-dev \  
                       libtinyxml-dev \  
                       qtbase5-dev \  
                       libvtk9-dev \  
                       libvtk9-qt-dev \  
                       gengetopt \  
                       help2man \  
                       groff \  
                       pod2pdf \  
                       bison \  
                       flex \  
                       libhpdf-dev \  
                       libtool  \
```

**Set Up Python Virtual Environment for OpenEMS**
Run the setup script and activate the virtual environment:
```bash
./setup_venv.sh.sh
```
 
Note the second time you can just run
```bash
source .venv/bin/activate
```

