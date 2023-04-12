# Some dummy code as we don't need to do anything to build simplemd 
# if plumed has been built
echo Hello simplemd world

# Check the required version of plumed exists
if ! command -v plumed$suffix &> /dev/null ; then
   echo install_plumed$suffix: broken >> $basedir/tests/simplemd/info.yml
else
   echo install_plumed$suffix: working >> $basedir/tests/simplemd/info.yml
fi
