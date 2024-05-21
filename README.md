# Complementarity

This is a jupyter notebook and associated python module to explore quantum complementarity from a discrete signal processing point of view. `Complementarity' is another way to refer to the way in which the uncertainty principle applies to physical observables that satisfy a certain non-commutative algebra. Here, we are using the interpretation from quantum mechanics where we think of observable quantities as Hermitian operators that act on a space of wavefunctions. A particular wavefunction has different representations but is in its essence the same object (a direction in a Hilbert space).

This particular notebook approaches complementarity from the engineering folk-wisdom that the time and frequency spreads of the representations of a particular signal show inverse proportionalities (roughly speaking). This project tries to study this trade-off quantitatively using discrete signal processing. The representation consists of the weights at each discrete sample location. It is intended to be useful both for teaching and exploration for interested folks. 

Files:

wavefunction.py - python module containing eponymous class for treating a particular 1-dimensional discretized setup using position as well as momentum representations (and associated operators). It allows a few different function families to be studied for the complementarity of their representations. A sweep function is added to study parameter-dependent variation.

complementarity.ipynb - Jupyter notebook which uses the above module to investigate different function families for the uncertainty principle and obtain functions that saturate the lower bound.

The term `complementarity' is based on Niels Bohr's article in Philosophy of Science from 1937 (and later in Science in the 50's):

"I am afraid that the short indications to which I have been obliged to restrict myself with respect to the last and many other points of this lecture will remind you only too well that in the last resort the direct use of any word must stand in complementary relationship to an analysis of its meaning. I hope, however, that I have to some extent succeeded in giving you the impression that my attitude is in no way in conflict with our common endeavors to arrive at as great a unification of knowledge as possible by the combating of prejudices in every field of research."

## Notation and key results

$$ | \psi \rangle = \Sigma_x \psi_x | x \rangle $$

![Wavefunctions illustrating the uncertainty tradeoff in position and momentum representations](20240521_illustration.png)



## References 

[1] [N. Bohr, "Causality and complementarity", Philosophy of Science, vol. 4, no. 3, July 1937](https://www.cambridge.org/core/journals/philosophy-of-science/article/abs/causality-and-complementarity/C193DEAB5C18330DD3739664761E8ECE)

[2] [N. Bohr, "On the notions of causality and complementarity", Science, Vol. 111, January 1950](https://www.science.org/doi/abs/10.1126/science.111.2873.51)

[3] [R. Shankar, "Principles of quantum mechanics", 2nd edition, Chapter 9: the Heisenberg Uncertainty Relations (pg. 237-241)]

[4] [M. A. Nielsen and I. C. Chuang, "Quantum computation and quantum information", Chapter 2: The postulates of quantum mechanics (pg. 89)]

[5] [Related work: S. Massar and P. Spindel, "Uncertainty relation for the discrete Fourier transform", Phys. Rev. Lett., 100, 190401, 2008](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.190401)

