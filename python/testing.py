import params;
import sims;

sims.dyeRaytracingNoClad("Rh6G", 1.5e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("Rh6G", 3e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("RhB", 3e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("LumR", 3e22, 1e-3, 0.97, 0.1);
sims.dyeRaytracingInClad("LumO", 3e22, 1e-3, 0.97, 0.1);
