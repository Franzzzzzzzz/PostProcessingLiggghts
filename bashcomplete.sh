_PPcomplete()
{
  local cur prev opts
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"
  opts="--setRadius --setRhograin --setG --setD --setR --setL --setLSimu --setDSimu --writing --W --alias --vtk --separe --w/forcetot --w/coupletot --w/otige --group-triangles --anglebystep --recode --norebuild --cp --uncompress --coarse-graining --coarse-graining-basic --substract --delcyl --deltank --set-vy-zero --use-box --setdeltaboite --sig-boites --winboxyper --wingauss --wingauss3D --w/kinetic --pressureFgrain --surface --surface2D --surfaces --xray --dumpall --downsampling --dstminmax --extract --mean --wallforce-by-angle --dump2restart --dump2vtk --w/speed --w/force --w/id --w/rayon --w/masse --dump2vtk --compress --reechantillonner --chainforce --donotusecfpos --forcetank --filter --F --filterCF --FCF --wallchainforce --wallchainforcenodump --clone-periodic-chain --cutoff --noperiodicity --v1 --grainforce-by-angle --grainforce --grain-rayon-around --is2D --break --noatmcoarse --wincreneau --iscylperiodic --mkdir --justtmp --forcetmp --nofileerror --clean --restart2vtk --help --bashcomplete --version --teststress "


  case "${prev}" in
    --setRadius) 
      COMPREPLY=( [setRadius] valeur )
      return 0
      ;;
    --setRhograin) 
      COMPREPLY=( [setRhograin] valeur )
      return 0
      ;;
    --setG) 
      COMPREPLY=( [setG] valeur )
      return 0
      ;;
    --setD) 
      COMPREPLY=( [setD] valeur )
      return 0
      ;;
    --setR) 
      COMPREPLY=( [setR] valeur )
      return 0
      ;;
    --setL) 
      COMPREPLY=( [setL] valeur )
      return 0
      ;;
    --setLSimu) 
      COMPREPLY=( [setLSimu] valeur )
      return 0
      ;;
    --setDSimu) 
      COMPREPLY=( [setDSimu] valeur )
      return 0
      ;;
    --anglebystep) 
      COMPREPLY=( [anglebystep] anglebystep )
      return 0
      ;;
    --coarse-graining) 
      COMPREPLY=( [coarse-graining] nbbox_x nbbox_y nbbox_z )
      return 0
      ;;
    --coarse-graining-basic) 
      COMPREPLY=( [coarse-graining-basic] nbbox_x nbbox_y nbbox_z )
      return 0
      ;;
    --substract) 
      COMPREPLY=( [substract] soustraire )
      return 0
      ;;
    --delcyl) 
      COMPREPLY=( [delcyl] xcyl ycyl zcyl dirxcyl dirycyl dirzcyl rayon longueur )
      return 0
      ;;
    --deltank) 
      COMPREPLY=( [deltank] rayon )
      return 0
      ;;
    --use-box) 
      COMPREPLY=( [use-box] box_xmin box_xmax box_ymin box_ymax box_zmin box_zmax )
      return 0
      ;;
    --setdeltaboite) 
      COMPREPLY=( [setdeltaboite] valeur )
      return 0
      ;;
    --sig-boites) 
      COMPREPLY=( [sig-boites] sig-x sig-y sig-z )
      return 0
      ;;
    --wingauss) 
      COMPREPLY=( [wingauss] sigma )
      return 0
      ;;
    --wingauss3D) 
      COMPREPLY=( [wingauss3D] sigmax sigmay sigmaz )
      return 0
      ;;
    --surface) 
      COMPREPLY=( [surface] nbbox_x nbbox_y )
      return 0
      ;;
    --surface2D) 
      COMPREPLY=( [surface2D] nbbox_x )
      return 0
      ;;
    --surfaces) 
      COMPREPLY=( [surfaces] nbbox_x nbbox_y )
      return 0
      ;;
    --xray) 
      COMPREPLY=( [xray] dir width height )
      return 0
      ;;
    --downsampling) 
      COMPREPLY=( [downsampling] downsampling )
      return 0
      ;;
    --dstminmax) 
      COMPREPLY=( [dstminmax] timestep )
      return 0
      ;;
    --extract) 
      COMPREPLY=( [extract] extract_deb extract_fin )
      return 0
      ;;
    --wallforce-by-angle) 
      COMPREPLY=( [wallforce-by-angle] xcyl zcyl nbbox_theta sigma )
      return 0
      ;;
    --dump2restart) 
      COMPREPLY=( [dump2restart] timestep )
      return 0
      ;;
    --forcetank) 
      COMPREPLY=( [forcetank] slices )
      return 0
      ;;
    --wallchainforce) 
      COMPREPLY=( [wallchainforce] xcyl zcyl rayon )
      return 0
      ;;
    --cutoff) 
      COMPREPLY=( [cutoff] cutoff )
      return 0
      ;;
    --grainforce-by-angle) 
      COMPREPLY=( [grainforce-by-angle] id nbbox_theta sigma )
      return 0
      ;;
    --grainforce) 
      COMPREPLY=( [grainforce] id )
      return 0
      ;;
    --grain-rayon-around) 
      COMPREPLY=( [grain-rayon-around] id )
      return 0
      ;;
    *)
      ;;
  esac

  if [[ ${cur} == -* ]] ; then
  COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
  return 0
  fi

  COMPREPLY=( $(compgen -f ${cur}) ) 
}