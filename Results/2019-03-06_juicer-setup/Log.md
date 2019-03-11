# 2019-03-06

First attempt at installing juicer (I feel like this is going to take a few tries).

```shell
# clone repo
cd $J
git clone git@github.com:theaidenlab/juicer.git

# go to project directory
cd $J/Davos/Src/
mkdir Juicer
cd Juicer

# make directory structure
ln -s $J/juicer/CPU scripts

# setup juicer_tools
pushd scripts/common
curl -O http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.8.9_jcuda.0.8.jar
ln -s juicer_tools.1.8.9_jcuda.0.8.jar juicer_tools.jar
popd

# setup references directory
mkdir references
```

This is getting crazy and hard to follow.
I'm not interested in continuing down this rabbit hole, I'm going to look elsewhere for analysis tools.
