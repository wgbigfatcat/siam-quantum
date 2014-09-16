QMD - Quantum Moleculear Dynamics is a program written by Chutchawan Jaisuk
during his undergraduate senior project in 2011. I have included this source
code as a starting point for further development in the future.

In his undergraduate senior project, Chutchawan studied a system of 28
water molecules using quantum molecular dynamics simulation. At time=0, the 
position of every atom in the system was initialized using "solven" module 
in VMD. Since, this was only the initial guess, no rigorous treatment was 
needed. Then, for every time step of the molecular dynamics simulation, 
a quantum mechanical calculation using Hartree-Fock method was done in order
to compute the net force acting on each nucleus. Once the net forces were found,
the next position of the nuclei in the next time step can be computed using
simple Newtonian mechanics. This procedure was sometimes called
Born-Oppenheimer Quantum Molecular Dynamics (as supposed to Car-Parrinello 
Quantum Molecule Dynamics).

After 1,000 time step, he was able to compute a certain statistical properties
of the system such as oxygen-oxygen radial distribution, which I was pround to
say, was in good agreement with the experimental distribution.

What make his senior project unique is how he handles the large number of 
water molecules in the system. Instead of calculating the entire system at once 
quantum mechanically to compute the net force acting on each nuecleus, he 
divides the 28 molecules into pairs, and computes interaction between pairs.
With careful book keeping, he is able to construct the net force from the 
force beteen pairs. In this way, he avoids spending a lot of time compute 
quantum mechanically on large system. This method is sometimes called pairwise
approximation, and it proves to be quite a good approximation for studying
water system.

You can read more about this topic and his works at the following website:

http://www.physics.kku.ac.th/Computational_physics (Selected Publications)
http://www.physics.kku.ac.th/sq

But be warned, his report is in Thai language. Feel free to send me email to 
ask questions about the project.

Chutchawan Jaisuk is now graduated and, at the time I am writing this report, 
he is working on Ph.D. degree in Mahidol University in Thailand.

Teepanis Chachiyo, June 2011. (teepanis@kku.ac.th) / (teepanis@hotmail.com)

