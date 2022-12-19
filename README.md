# ringSNARK - A C++ Library for zkSNARK arguments over ring
This repository contains an implementation of the Rinocchio [[1]](#1) protocol (SNARK for Ring Arithmetic) for general rings. 
This repository also contains instantiates the general protocol for rings $Z_{\prod_i q_i}[X]/\langle X^N+1\rangle$ as used in Fully Homomorphic Encryption (FHE), notably those used by the SEAL [[2]](#2) library. 


## About
ringSNARK is being developed in the [Privacy-Preserving Systems Lab](https://pps-lab.com) at [ETH Zurich](https://ethz.ch/en.html) by [Christian Knabenhans](https://cknabs.github.io). 

For its concrete instantiation for FHE rings, ringSNARK uses the [Microsoft SEAL](https://github.com/microsoft/SEAL) library and the `polytools` set of utilities. `polytools` was developed by [Alexander Viand](https://pps-lab.com/people/alexanderviand/) at Intel Labs, and is still closed-sourced (stay tuned for its open-sourcing!). 

### Rinocchio: SNARKs for Ring Arithmetic
This code implements the "Rinocchio" protocol outlined in [[1]](#1); [Rinocchio.pdf](/blob/master/Rinocchio.pdf) contains an implementation-friendly description of the protocol, a concrete runtime analysis, and our optimizations to the protocol (in particular, faster encodings for FHE rings). 

## Security

The theoretical security of the underlying SNARKs and their assumptions are analyzed in [[1]](#1).  
This code is a research-quality proof-of-concept, has not undergone a thorough security review, and is still being actively developed. 
You are welcome to use it for proof-of-concept and academic projects, but this code is not suitable for critical and production systems. 

## Build instructions

The ringSNARK library relies on the following:
- C++ build environment
- CMake build infrastructure

Optionally, if support for FHE rings is needed, the following are required: 
- SEAL
- polytools

### Building
```bash
git clone https://github.com/MarbleHE/ringSNARK && cd ringSNARK
git submodule init && git submodule update --recursive
mkdir build && cd build && cmake ..
make
```

## Directory structure 
`include/` contains header and template files for the Rinocchio protocol over general rings, as well as header files for an instantiation of Rinocchio for FHE schemes as implemented in SEAL. 

`src/` contains the implementation of the SEAL-specific prover and verifier. 

`example.cpp` shows how to instantiate the code to prove and verify a small (1 addition, 2 multiplications) circuit over FHE rings. 

## Roadmap
This code is being actively developed, with major changes in the interfaces coming soon. Here is the short-term implementation plan: 

- [ ] Migrate to a libsnark-inspired interface
  - [ ] Implement circuit -> R1CS -> QRP generation pipeline
  - [ ] Change interfaces to expose same API as libsnark
  - [ ] Abstract random sampling inside keypair generation process
- [ ] Implement optimized FHE-ring encoding/decoding
- [ ] Add full zero-knowledge support
- [ ] Analyze and implement RingGroth16 in addition to Rinocchio

## References

<a id="1">[1]</a> C. Ganesh, A. Nitulescu, and E. Soria-Vazquez, Rinocchio: SNARKs for Ring Arithmetic. Cryptology ePrint Archive, Paper 2021/322, 2021. [Online]. Available: https://eprint.iacr.org/2021/322

<a id="2">[2]</a> Microsoft SEAL (release 4.0). https://github.com/Microsoft/SEAL, 2022. [Online]. Available: https://github.com/Microsoft/SEAL 
