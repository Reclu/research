#!/bin/sh

for FILE in ./*.csv ; do 
  cp "${FILE}" "${FILE::-4}.csv.org"
  sed 's/,/ /g' "${FILE}" >"${FILE}2"
  sed 's/"/ /g' "${FILE}2" >"${FILE}3"
  rm "${FILE}2"  "${FILE}"
  mv "${FILE}3" "${FILE}"
done
