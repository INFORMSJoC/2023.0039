[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An Optimization-based Scheduling Methodology for Appointment Systems with Heterogeneous Customers and Non-stationary Arrival Processes

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The data in this repository exemplifies the data
that was used in the research reported on in the paper 
[An Optimization-based Scheduling Methodology for Appointment Systems with Heterogeneous Customers and Non-stationary Arrival Processes](https://doi.org/10.1287/ijoc.2023.0039) by S. Chatterjee, Y. Hebaish, H. Aprahamian and L. Ntaimo. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0039

https://doi.org/10.1287/ijoc.2023.0039.cd

Below is the BibTex for citing this snapshot of the repository.

```
@article{chatterjee2024,
  author =        {Chatterjee, S. and Hebaish, Y. and Aprahamian, H. and Ntaimo, L.},
  publisher =     {INFORMS Journal on Computing},
  title =         {{An Optimization-based Scheduling Methodology for Appointment Systems with Heterogeneous Customers and Non-stationary Arrival Processes}},
  year =          {2024},
  doi =           {10.1287/ijoc.2023.0039.cd},
  url =           {https://github.com/INFORMSJoC/2023.0039},
  note =          {Available for download at https://github.com/INFORMSJoC/2023.0039},
}  
```

## Description

The [data folder](data) includes the arrival rate data for first-time and crisis customers and the optimal schedule proportions for both service types for all the cases in Table 1 and Figure 9 of the paper. It also contains the schedule proportions obtained using the DPS policy, as shown in Figure 10 of the paper. Refer to the [README.md file in that folder](data/README.md) for more details on the data.

The [scripts folder](scripts) includes all the Python codes required to generate the various scheduling polcies used in this paper, as well as the results in the Online Supplement. It also includes a 'requirements.txt' file that can be used to load the relevant packages needed to run the scripts into a virtual Python environment. Note that to run some of the scripts, a license of the Gurobi solver will be required. Refer to the [README.md file in that folder](scripts/README.md) for more details on the parameters needed for the scripts.

## Contact
If you have any question please contact me at:
- Sohom Chatterjee (sohom070994@gmail.com)
