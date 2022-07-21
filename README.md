# rinocchio-SEAL
This repository contains an implementation of the Rinocchio [[1]](#1) protocol (SNARK for Ring Arithmetic), in particular over SEAL[[2]](#2) FHE ring elements. 


## Structure 
`include/` contains header and template files for the Rinocchio protocol over general rings, as well as header files for an instantiation of Rinocchio for FHE schemes as implemented in SEAL. 

`src/` contains the implementation of the SEAL-specific prover and verifier. 

`example.cpp` shows how to instantiate the code to prove and verify a small (1 addition, 2 multiplications) circuit over FHE rings. 

## Roadmap

- [ ] Implement QRP generation from circuit/R1CS
- [ ] Add random sampling interface for keying material
- [ ] Implement encoding/decoding
- [ ] Add zero-knowledge

## License

TBD

## Contact

Christian Knabenhans - christian.knabenhans@alumni.ethz.ch

## References

<a id="1">[1]</a> C. Ganesh, A. Nitulescu, and E. Soria-Vazquez, Rinocchio: SNARKs for Ring Arithmetic. Cryptology ePrint Archive, Paper 2021/322, 2021. [Online]. Available: https://eprint.iacr.org/2021/322

<a id="2">[2]</a> Microsoft SEAL (release 4.0). https://github.com/Microsoft/SEAL, 2022. [Online]. Available: https://github.com/Microsoft/SEAL 
