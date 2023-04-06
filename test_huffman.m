%Cleans environment
clc;
clear;
close all;

%Number of symbols to simulate
N=1024;

%Number of bits
NB=11;

%Number of symbols
NS=2^NB;

%Zero probability
pz=0.7;

%Symbols probabilities
probs=[pz (1-pz)/(NS-1)*ones(1,NS-1)];

%Random symbols 
symbols = randsample([0:(NS-1)], NS, true, probs); 
 
%Gets huffman dictionary
[dict] = huffman_code([]  , probs);

%Huffman encoding
[encoded] = huffman_code(symbols,[], dict);

%Hufman decoding
decoded = huffman_decode(encoded, dict);

%l1-error
fprintf('l1-error: %f\n',sum(symbols-decoded));