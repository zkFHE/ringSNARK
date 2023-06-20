# ringSNARK - A modular library for zkSNARKs over rings
This repository contains a generic implementation of the Rinocchio and ringGroth16 [[1]](#1) SNARKs for general rings. 

This repository also contains concrete instantiations of encodings for the rings $Z_q$, $Z_q^N$, and $Z_q[X]/\langle X^N+1\rangle$ for a composite $q$. 
These rings are especially useful to (efficiently) prove statements about lattice relations, in particular for Fully Homomorphic Encryption (FHE). 

## About
### Frontend
This library provides a [libsnark](https://github.com/scipr-lab/libsnark)-inspired domain-specific language to build gadgets and specify constraints. 

### Proof systems
ringSNARK implements the two proof systems from [[1]](#1): Rinocchio (based on Pinocchio SNARK) and ringGroth16 (based on Groth16). 

For the ring $Z_q^N$, we use _batched_ encodings, which are orders of magnitude much more efficient than the ones proposed in [[1]](#1). 

### Backend
ringSNARK can use two backends for fast vector/polynomial ring arithmetic: 
- [Microsoft SEAL](https://github.com/microsoft/SEAL), via the [SEAL-Polytools](https://github.com/MarbleHE/SEAL-Polytools) arithmetic wrapper
- [OpenFHE](https://github.com/openfheorg/openfhe-development)
### Structure 
```
├ docs --------------- auxiliary material, including specifications, scripts, and presentations
├ examples ----------- circuits for various (FHE) use cases
└ ringsnark
  ├ gadgetlib -------- libsnark-style gadgets
  ├ reductions ------- libsnark-style for R1CS->QRP translation
  ├ relations -------- libsnark-style data structures for R1CS/QRP instances
  ├ zk_proof_systems - template implementation of Rinocchio and ringGroth16
  └ seal ------------- rings implemented with the SEAL backend
```

### Security
The theoretical security of the underlying SNARKs and their assumptions are analyzed in [[1]](#1).  
This code is a research-quality proof-of-concept, has not undergone a thorough security review, and is still being actively developed. 
You are welcome to use it for proof-of-concept and academic projects, but this code is not suitable for critical and production systems. 

## Build instructions

The ringSNARK library relies on the following:
- C++ build environment
- CMake build infrastructure

### Requirements
This library requires the `boost` C++ library. 

Optionally, if support for FHE rings is needed, the following dependencies are needed: 
- For the SEAL backend: [SEAL](https://github.com/microsoft/SEAL) (tested with versions 4.0.0 to 4.1.1) and [SEAL-Polytools](https://MarbleHE/SEAL-Polytools)
- For the OpenFHE backend: [OpenFHE](https://github.com/openfheorg/openfhe-development)
  
Both are fetched automatically as submodules, so you don't need to install them separately. 

### Building
```bash
git clone https://github.com/MarbleHE/ringSNARK && cd ringSNARK
git submodule init && git submodule update --recursive
mkdir build && cd build && cmake ..
make
```

## References
<a id="1">[1]</a> C. Ganesh, A. Nitulescu, and E. Soria-Vazquez, Rinocchio: SNARKs for Ring Arithmetic. Cryptology ePrint Archive, Paper 2021/322, 2021. Available: https://eprint.iacr.org/2021/322
