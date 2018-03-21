-- /* Copyright (c) 2016 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- DEAR LUA USER:
-- This is the developer version of the lua parameter file. 
-- It gives all the options you can have. 
-- Many of these the client will not need.

-- NOTE --
-- if you are using single component plummer model, it will take the baryonic
-- matter component parameters. meaning you input should look like
-- ft, bt, rscale_baryon, radius_ratio, baryon mass, mass ratio
-- typical parameters: 4.0, 1.0, 0.2, 0.2, 12, 0.2
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- STANDARD  SETTINGS   -- -- -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
totalBodies           = 20000   -- -- NUMBER OF BODIES           -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.68"  -- -- MINIMUM APP VERSION        -- --

run_null_potential    = false   -- -- NULL POTENTIAL SWITCH      -- --
two_component_model   = true    -- -- TWO COMPONENTS SWITCH      -- --
use_tree_code         = true    -- -- USE TREE CODE NOT EXACT    -- --
print_reverse_orbit   = false   -- -- PRINT REVERSE ORBIT SWITCH -- --
print_out_parameters  = false    -- -- PRINT OUT ALL PARAMETERS   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- PARAMETER SETTINGS   -- -- -- -- -- -- -- -- -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

-- -- -- -- -- -- -- -- -- HISTOGRAM   -- -- -- -- -- -- -- -- -- -- -- -- --
lda_bins        = 26      -- number of bins in lamdba direction
-- lda_lower_range = -150    -- lower range for lambda
-- lda_upper_range = 150     -- upepr range for lamdba
lda_lower_range = -38.235294117647058    -- lower range for lambda
lda_upper_range = 38.235294117647058     -- upepr range for lamdba

bta_bins        = 1       -- number of beta bins. normally use 1 for 1D hist
bta_lower_range = -15     -- lower range for beta
bta_upper_range = 15      -- upper range for beta

SigmaCutoff          = 2.5     -- -- sigma cutoff for outlier rejection -- --
Correction           = 1.111   -- -- correction for outlier rejection -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 

-- -- -- -- -- -- -- -- -- AlGORITHM OPTIONS -- -- -- -- -- -- -- --
use_best_likelihood  = true    -- use the best likelihood return code
best_like_start      = 0.98    -- what percent of sim to start

use_beta_disps       = true    -- use beta dispersions in likelihood
use_vel_disps        = false   -- use velocity dispersions in likelihood

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- ADVANCED DEVELOPER OPTIONS -- -- -- -- -- -- -- --        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- These options only work if you compile nbody with  -- -- --
-- -- -- -- -- -- the -DNBODY_DEV_OPTIONS set to on                  -- -- --   

useMultiOutputs       = false   -- -- WRITE MULTIPLE OUTPUTS       -- --
freqOfOutputs         = 6       -- -- FREQUENCY OF WRITING OUTPUTS -- --

timestep_control     = false    -- -- control number of steps      -- --
Ntime_steps          = 1        -- -- number of timesteps to run   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        




-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
orbit_parameter_l  = 218
orbit_parameter_b  = 53.5
orbit_parameter_r  = 28.6
orbit_parameter_vx = -156 
orbit_parameter_vy = 79 
orbit_parameter_vz = 107
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        

        
function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        return  Potential.create{
            spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
            disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
            halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
        }
   end
end


function get_timestep()
    if(two_component_model == true) then
        --Mass of a single dark matter sphere enclosed within light rscale
        mass_enc_d = mass_d * (rscale_l)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

        --Mass of a single light matter sphere enclosed within dark rscale
        mass_enc_l = mass_l * (rscale_d)^3 * ( (rscale_l)^2 + (rscale_d)^2  )^(-3.0/2.0)

        s1 = (rscale_l)^3 / (mass_enc_d + mass_l)
        s2 = (rscale_d)^3 / (mass_enc_l + mass_d)
        
        --return the smaller time step
        if(s1 < s2) then
            s = s1
        else
            s = s2
        end
        
        -- I did it this way so there was only one place to change the time step. 
        t = (1 / 100.0) * ( pi_4_3 * s)^(1.0/2.0)
        
    --     tmp = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_d)) / (mass_l + mass_d))
    --     print('timestep ', t, tmp)
    else 
        t = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_l)) / (mass_l))
    end
--     print("LUA FILE: ", t)
    return t
end


function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timestep    = get_timestep(),
      eps2        = calculateEps2(totalBodies, soften_length ),
      criterion   = criterion,
      useQuad     = true,
      useBestLike   = use_best_likelihood,
      BestLikeStart = best_like_start,
      useBetaDisp   = use_beta_disps,
      useVelDisp    = use_vel_disps,
      Nstep_control = timestep_control,
      Ntsteps       = Ntime_steps,
      MultiOutput   = useMultiOutputs,
      OutputFreq    = freqOfOutputs,
      BetaSigma     = SigmaCutoff,
      VelSigma      = SigmaCutoff,
      BetaCorrect   = Correction,
      VelCorrect    = Correction,
      theta       = 1.0
   }
end


function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity, finalPosition2, finalVelocity2
    if(run_null_potential == true) then
        print("placing dwarf at origin")
        finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
    else 
        finalPosition, finalVelocity = reverseOrbit{
            potential = potential,
            position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
            velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
            tstop     = revOrbTime,
            dt        = ctx.timestep / 10.0
            }
            
--         finalPosition2, finalVelocity2 = reverseOrbit2{
--             potential = potential,
--             position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
--             velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
--             tstop     = revOrbTime,
--             dt        = ctx.timestep / 10.0
--             }
--         print(finalPosition, finalVelocity)
--         print(finalPosition2, finalVelocity2)
    end
    
    if(print_reverse_orbit == true) then
        local placeholderPos, placeholderVel = PrintReverseOrbit{
            potential = potential,
            position  = lbrToCartesian(ctx, Vector.create(orbit_parameter_l, orbit_parameter_b, orbit_parameter_r)),
            velocity  = Vector.create(orbit_parameter_vx, orbit_parameter_vy, orbit_parameter_vz),
            tstop     = .14,
            tstopf    = .20,
            dt        = ctx.timestep / 10.0
        }
        print('Printing reverse orbit')
    end
    

  
    if(two_component_model) then 
        firstModel = predefinedModels.mixeddwarf{
            nbody       = totalBodies,
            prng        = prng,
            position    = finalPosition,
            velocity    = finalVelocity,
            comp1       = Dwarf.plummer{mass = mass_l, scaleLength = rscale_l}, -- Dwarf Options: plummer, nfw, general_hernquist
            comp2       = Dwarf.plummer{mass = mass_d, scaleLength = rscale_d}, -- Dwarf Options: plummer, nfw, general_hernquist
            ignore      = true
        }
       
    else
        firstModel = predefinedModels.plummer{
            nbody       = totalBodies,
            prng        = prng,
            position    = finalPosition,
            velocity    = finalVelocity,
            mass        = mass_l,
            scaleRadius = rscale_l,
            ignore      = false
        }
  
    end
  
  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     
     -- ANGULAR RANGE AND NUMBER OF BINS
     lambdaStart = lda_lower_range,
     lambdaEnd   = lda_upper_range,
     lambdaBins  = lda_bins,
     
     betaStart = bta_lower_range,
     betaEnd   = bta_upper_range,
     betaBins  = bta_bins
}
end


arg = { ... } -- -- TAKING USER INPUT
assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
-- argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
argSeed = 7854614814 -- -- SETTING SEED TO FIXED VALUE Seed used for mw@home
prng = DSFMT.create(argSeed)

-- -- -- -- -- -- -- -- -- ROUNDING USER INPUT -- -- -- -- -- -- -- --
function round(num, places)
  local mult = 10.0^(places)
  return floor(num * mult + 0.5) / mult
end

-- -- -- -- -- -- ROUNDING TO AVOID DIFFERENT COMPUTER TERMINAL PRECISION -- -- -- -- -- --
dec = 9.0
evolveTime       = round( tonumber(arg[1]), dec )
rev_ratio        = round( tonumber(arg[2]), dec )
rscale_l         = round( tonumber(arg[3]), dec )
light_r_ratio    = round( tonumber(arg[4]), dec )
mass_l           = round( tonumber(arg[5]), dec )
light_mass_ratio = round( tonumber(arg[6]), dec )

-- -- -- -- -- -- -- -- -- DWARF PARAMETERS   -- -- -- -- -- -- -- --
revOrbTime = evolveTime
dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

if(use_tree_code) then
    criterion = "TreeCode"
--     criterion = "NewCriterion"
else
    criterion = "Exact"
end

if(print_out_parameters) then
    print('forward time=', evolveTime, '\nrev time=',  revOrbTime)
    print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
    print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
    print('total mass solar= ', (mass_d + mass_l) * 222288.47)
    print('rl = ', rscale_l, 'rd = ', rscale_d)
end