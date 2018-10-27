#!/usr/bin/env bash

# Credenciales para conectarse al servidor en el que se despliega el contenedor;


# Debe estar definido en Gitlab, como una variable del proyecto


DEPLOYMENT_PASS=$1

IP=quantum@10.0.0.240

# Conectarse por SSH al servidor de builds, 

sshpass -p "$DEPLOYMENT_PASS" ssh -o StrictHostKeyChecking=no -T "$IP" << EOF
mkdir -p /home/quantum/QSimulator
cd /home/quantum/QSimulator
rm -rf *
git init
git remote add origin http://10.0.0.250/quantum/QSimulator.git
git pull -u origin master
exit
EOF