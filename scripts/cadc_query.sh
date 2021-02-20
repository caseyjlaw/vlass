#!/bin/bash
RA=$1
DEC=$2
RADIUS=$3

CUTOUT_SERVICE="http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/cutout"

QUERY="SELECT TOP 10 Artifact.uri
FROM caom2.Observation AS Observation 
   JOIN caom2.Plane AS Plane ON Observation.obsID=Plane.obsID 
   JOIN caom2.Artifact AS Artifact ON Artifact.planeID=Plane.planeID
WHERE collection='VLASS' 
   AND ( INTERSECTS( CIRCLE('ICRS', $RA, $DEC, $RADIUS), Plane.position_bounds ) = 1 ) 
   AND Artifact.productType='science'"

for uri in `cadc-tap query -s ivo://cadc.nrc.ca/tap -f tsv "$QUERY"`;
do 
   echo ${uri} | grep "ad"  || continue
   if [ ! -f `echo ${uri} | grep "ad" | cut -d '/' -f 2` ];
   then
      uri=`echo ${uri} | sed -e 's/+/%2B/'`
      data="uri=${uri}&cutout=Circle+ICRS+${RA}+${DEC}+${RADIUS}"
      curl -L -O  ${CUTOUT_SERVICE}?${data}
   fi
done
