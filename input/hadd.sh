#!/bin/sh

# Arguments and ch setting
filename=$1
ch=`echo $filename | cut -d"." -f4`

# Get sample list
data=(` ls $filename | grep data | sed -e 's/submitDir_//g'`)
Wjets=(` ls $filename | grep W | grep Sh221$ | grep -v Wqq | grep -v Wlv | sed -e 's/submitDir_//g'`)
Zjets=(` ls $filename | grep Z | grep Sh221$ | grep -v Wqq | grep -v Zqq | sed -e 's/submitDir_//g'`)
ttbar=(` ls $filename | grep -e ttbar -e sto | sed -e 's/submitDir_//g'`)
diboson=(` ls $filename | grep -e Wqq -e Zqq -e Wlv -e Zll -e Zvv | grep -v Hbb | sed -e 's/submitDir_//g'`)
if [ "$ch" = "0Lep" ] ; then
	tagW='HVT'
	tagZ='RS'
elif [ "$ch" = "1Lep" ] ; then
	tagW='HVTWW'
	tagZ='HVT'
elif [ "$ch" = "2Lep" ] ; then
	tagW='HVT'
	tagZ='RS'
fi
signalW=(` ls $filename | grep  $tagW | sed -e 's/submitDir_//g'`)
signalZ=(` ls $filename | grep  $tagZ | sed -e 's/submitDir_//g'`)
signal=('sigDBL')

# Directory for output
if test ! -d merged ; then
	mkdir merged merged/allMC merged/Wjets merged/Zjets merged/diboson merged/ttbar merged/hbb merged/signalW merged/signalZ merged/data merged/sigDBL
fi
rm -rf merged/allMC/*$ch*

#<<COMMENTOUT
# hadd for each sampl
echo '++ data ++++++++++++++++++++++++++++++'
rm -rf merged/data/*$ch*
for i in ${data[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/data/$i'.'$ch'.root'
done
cd merged/data
ls -S | grep $ch | xargs hadd data.$ch'.root'
ls | grep -v 'data.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../


#COMMENTOUT
echo '++ Wjets ++++++++++++++++++++++++++++++'
rm -rf merged/Wjets/*$ch*
for i in ${Wjets[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/Wjets/$i'.'$ch'.root'
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/allMC/$i'.'$ch'.root'
done
cd merged/Wjets
ls -S | grep $ch | xargs hadd Wjets.$ch'.root'
ls | grep -v 'Wjets.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../

#COMMENTOUT
echo '++ Zjets ++++++++++++++++++++++++++++++'
rm -rf merged/Zjets/*$ch*
for i in ${Zjets[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/Zjets/$i'.'$ch'.root'
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/allMC/$i'.'$ch'.root'
done
cd merged/Zjets
ls -S | grep $ch | xargs hadd Zjets.$ch'.root'
ls | grep -v 'Zjets.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../

#<<COMMENTOUT
echo '++ diboson ++++++++++++++++++++++++++++++'
rm -rf merged/diboson/*$ch*
for i in ${diboson[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/diboson/$i'.'$ch'.root'
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/allMC/$i'.'$ch'.root'
done
cd merged/diboson
ls -S | grep $ch | xargs hadd diboson.$ch'.root'
ls | grep -v 'diboson.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../

echo '++ ttbar ++++++++++++++++++++++++++++++'
rm -rf merged/ttbar/*$ch*
for i in ${ttbar[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/ttbar/$i'.'$ch'.root'
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/allMC/$i'.'$ch'.root'
done
cd merged/ttbar
ls -S | grep $ch | xargs hadd ttbar.$ch'.root'
ls | grep -v 'ttbar.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../

echo '++ signal W ++++++++++++++++++++++++++++++'
rm -rf merged/signalW/*$ch*
for i in ${signalW[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/signalW/$i'.'$ch'.root'
done
cd merged/signalW
ls -S | grep $ch | xargs hadd signalW.$ch'.root'
cd ../../

echo '++ signal Z ++++++++++++++++++++++++++++++'
rm -rf merged/signalZ/*$ch*
for i in ${signalZ[@]};do
	cp -r $filename/submitDir_$i/data-MVATree/*root merged/signalZ/$i'.'$ch'.root'
done
cd merged/signalZ
ls -S | grep $ch | xargs hadd signalZ.$ch'.root'
cd ../../

echo '++ ALL Bkg ++++++++++++++++++++++++++++++'
cd merged/allMC
ls -S | grep $ch | xargs hadd allMC.$ch'.root'
ls | grep -v 'allMC.'$ch'.root' | grep $ch | xargs rm -rf
cd ../../

