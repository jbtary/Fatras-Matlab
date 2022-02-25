clc
clear all
close all


%Llama la funcion GravFunc() con la palabra Model para llamar 
%los archivos de entrada Mdel_mod.txt y Model_X.txt
root='Model';
ax_ymin=-375e3;
ax_xmax=375e3;

%Llama a funcion que calcula la anomalia gravitacional
GravFunc(root,ax_ymin,ax_xmax)



