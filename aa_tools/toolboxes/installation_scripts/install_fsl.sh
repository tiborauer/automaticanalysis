# This is a script to download and install FSL. The script checks the installation by looking for the "fsl" folder
#
# source install_fsl.sh <installation_directory> <version> <os> <standard_only> <config_file>
#
# Input aguments
#   - installation_directory:   This is where freesurfer will be installed. When configuring aap.directory_conventions.fsldir MUST point to '<installation_directory>/fsl'
#   - version:                  Version number of the desired installation according to semantic versioning (x.y.z)
#   - standard_only:            Only fsl/data/standard is installed. The whole FSL is not installed.
#   - config_file:              Full path to the FSL configuration script to be generated

INSTDIR=$1
VERSION=$2
STANDARDONLY=$3
CONFIGFILE=$4

function fslconfig_bash {
    echo 'IFS=":"'
    echo 'read -ra pathvar <<< $PATH'
    echo 'condaPath=(); for p in ${pathvar[@]}; do if [[ $p == *"conda"* ]]; then condaPath+=($p); fi; done'
    echo 'export FSLDIR=INSTDIR/fsl'
    echo 'export PATH=$FSLDIR/bin:$PATH'
    echo '. $FSLDIR/etc/fslconf/fsl.sh'
    echo 'if [[ ${#condaPath[@]} > 0 ]]; then export PATH=$(printf "%s:" "${condaPath[@]}")$PATH; fi'
}

cd ${INSTDIR}
echo "Installing FSL ${VERSION} in ${INSTDIR}"

if [[ "x${STANDARDONLY}x" == "x1x" ]]; then
    wget -q http://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-6.0.5.2-centos7_64.tar.gz
    tar xzf fsl-6.0.5.2-centos7_64.tar.gz fsl/data/standard
    rm fsl-6.0.5.2-centos7_64.tar.gz
else
    wget -q https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/releases/fslinstaller.py
    python fslinstaller.py -V $VERSION -d ${INSTDIR}/fsl -s 1 -m 1

    # config script
    if [[ $(basename $(echo $SHELL)) == "bash" ]]; then
        CFGSTR=$(fslconfig_bash)
        echo "${CFGSTR//INSTDIR/$INSTDIR}" > $CONFIGFILE
    else
        echo "FSL configuration is implemented only for BASH."
        echo "FSL has not been installed properly. Exiting..." >&2
        exit -1; 
    fi

    rm ${INSTDIR}/fslinstaller.py
fi

if [[ -d "${INSTDIR}/fsl" ]]; then
    echo "FSL has been installed." >&2
else
    echo "FSL has not been installed properly. Exiting..." >&2
    exit -1; 
fi