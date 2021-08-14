# Codes-Janbazi-JPhysChemA
Supplemental Material for J Phys Chem A

                 
This repository is part of the Supplemental Material for
Thermochemistry of Oxygen-Containing Organosilane Radicals and Uncertainty Estimations of Organosilane Group-Additivity Values      
                                                                              
H. Janbazi 1, C. Schulz 2,3, I. Wlokas 1,3, S. Peukert 2,3                     
                                                                              
1 Institute of Combustion and Gas Dynamics (IVG)- Fluid Dynamics,                        
2 Institute of Combustion and Gas Dynamics (IVG)- Reactive Fluids,           
3 Center for Nanointegration Duisburg-Essen (CENIDE),                        
University of Duisburg-Essen, 47058 Duisburg, Germany                      
                                                                              
Journal of Physical Chemistry A, 2021                      


###################################################

To facilitate the application of GAVs of Si-organic species, GAVs derived from the present and previous works are summarized in an XML file and a Python44 script (Supplementary Material) was written that uses the GAVs stored in the XML file “groups_data” to extract ΔHf° and S° (in kJ/mol and J/(mol K)) as well as the NASA polynomials for thermochemical data in Chemkin format. Seven-coefficient NASA polynomials for the low temperature range from 298–1000 K and seven-coefficient polynomials for the high-temperature range from 1000–2000 K were obtained. The syntax for running the python script is: thermoCreator.py 'Molecule_name' 'GroupA:a GroupB:b GroupC:c ....'. In this syntax, a, b, and c represent the numbers how often the group increments A, B, and C appear in that molecule. The following example is given: (CH3)SiH2-(CH3)(SiH3)C•. This radical contains the C−(H)3(Si), Si−(RC)(C)(H)2, C−(RC)(H)3, Si−(RC)(H)3, and RC−(C)(Si)2 groups once per radical. In this case, the syntax is: thermoCreator.py '(CH3)SiH2-(CH3)(SiH3)C_radical' 'C-(H)3(Si):1 Si-(RC)(C)(H)2:1 C-(RC)(H)3:1 Si-(RC)(H)3:1 RC-(C)(Si)2:1'. The output for ΔHf° and S° is 104.2 kJ/mol and 421.7 J/(mol K). The python script and the XML file with GAVs are part of the supplementary information. 

                                                                              
                                                                              
                                                                              

