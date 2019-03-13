# Engine App Interface

Althogh we can get the symbolic representation of the QuantumLattices by our *unitcell-description* framework, there still remains a long way to implement concrete algorithms such as TBA, ED, etc. Despite the quite different techinal details, algorithms shares common functionalities to be furnished with:
* provide tasks to be conducted with controlling parameters;
* record the results of some tasks for later use or analysis;
* update some parameters of the QuantumLattices to reconduct tasks;
* keep logs during code executions for debug;
* cahe intermediate data to improve efficienty;
* ...
We thus provide a set of generic interfaces to resolve these problems, basiscally in the so called *Engine-App* mode. Specifically, algorithms are treated as `Engine` and taks as `App`. `Engine` deals with the cores of algorithms along with file mangament, parameter updating and data caching, while `App` decides the concrete tasks to be conducted and provides hyper controlling parameters.
