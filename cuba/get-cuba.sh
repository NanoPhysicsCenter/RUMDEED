#!/bin/bash
echo "Do you wish to download the CUBA library?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) wget https://feynarts.de/cuba/Cuba-4.2.2.tar.gz; break;;
        No ) echo "Please download CUBA from https://feynarts.de/cuba"; echo "Install it manually or place the tar.gz file in the cuba directory"; exit;;
    esac
done

