doMuonSel=1
doMigration=0
doEffPurity=0

Particle=NP
Neutrino="minervameNu"
playlists="minervame1A minervame1B minervame1C minervame1D minervame1E minervame1F minervame1G minervame1L minervame1M minervame1N minervame1O minervame1P"
multiplicity=2


source ../include/scripts_source_nu.sh
use2D=""

#for tune in pionrpaElasticNuwroSF CVElasticNuwroSF  CVElasticNuwroSFn CVElasticNuwroSFn2 CV2ElasticNuwroSF CV2ElasticNuwroSFn CV2ElasticNuwroSFn2 CV2Elastic CV2
for tune in  CV2Elastic
do
  fbase0="../rootfiles/${Neutrino}/cuts${cuts}${use2D}/${tune}${tag}/"
  mkdir -p $fbase0
  echo $fbase0
  fbaseorig=$basefolder$use2D/${Particle}_${tune}/
  ln -sf $basefolder${use2D}/${Particle}_${tune}/MuonEvent*root $fbase0
  #for sample in Signal BlobSideBand MichelSideBand MicBlobSideBand
  for sample in  Signal BlobSideBand MichelSideBand MicBlobSideBand
  do
    #fbase="$fbase0/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_"
    fbase="$fbase0/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_"
    cmd="./PlaylistAdder $fbase MuonEventSelection 0 $playlists"
    #for pl in $playlists
    #do
    #  fnamec="${fbaseorig}/MuonEventSelection_MakeFlux-1_Multiplicity-${multiplicity}_Sample-${sample}_$pl.root"
    #  #echo check $fnamec
    #  if test ! -f "$fnamec"
    #  then
    #    echo $pl does not exists
    #  fi
    #done 
    if [ $doMuonSel = 1 ]
    then
      echo $cmd
      nohup $cmd > nohup_nu_$sample.log 2>&1 &
    fi
  done


  mkdir -p $fbase0
  echo $fbase0
  for pl in $playlists
  do
    ln -sf $basefolder/${Particle}_${tune}/Migration_*${pl}_Q2.root $fbase0/Migration_${pl}.root
    ln -sf $basefolder/${Particle}_${tune}/EffPurity*${pl}*.root $fbase0/EffPurity_${pl}.root
  done

  fbase="$fbase0/EffPurity_"
  cmd="./PlaylistAdder $fbase EffPurity 0 $playlists"
  if [ $doEffPurity = 1 ]
  then
    echo $cmd
    $cmd 
  fi


  fbase="$fbase0/Migration_"
  cmd="./PlaylistAdder $fbase Migration 0 $playlists"
  if [ $doMigration = 1 ]
  then
    echo $cmd
    $cmd 
  fi

done


