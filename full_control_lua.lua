-- /* Copyright (c) 2016 Siddhartha Shelton */

-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- DEAR LUA USER:
-- This is the developer version of the lua parameter file. 
-- It gives all the options you can have. 
-- Many of these the client will not need.
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
        
        
        
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
totalBodies           = 20000   -- -- NUMBER OF BODIES           -- --
nbodyLikelihoodMethod = "EMD"   -- -- HIST COMPARE METHOD        -- --
nbodyMinVersion       = "1.64"  -- -- MINIMUM APP VERSION        -- --

run_null_potential    = false   -- -- NULL POTENTIAL SWITCH      -- --
two_component_model   = true    -- -- TWO COMPONENTS SWITCH      -- --
print_reverse_orbit   = false   -- -- PRINT REVERSE ORBIT SWITCH -- --
print_out_parameters  = true    -- -- PRINT OUT ALL PARAMETERS   -- --
-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- 
-- -- -- -- -- -- -- -- -- DWARF STARTING LOCATION   -- -- -- -- -- -- -- --
l  = 218
b  = 53.5
r  = 28.6
vx = -156 
vy = 79 
vz = 107
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
    
    tmp = sqr(1/10.0) * sqrt((pi_4_3 * cube(rscale_d)) / (mass_l + mass_d))
--     print('timestep ', t, tmp)
    
    return t
end


function makeContext()
   soften_length  = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve  = evolveTime,
      timestep    = get_timestep(),
      eps2        = calculateEps2(totalBodies, soften_length ),
      criterion   = "NewCriterion",
      useQuad     = true,
      useBestLike = true,
      useVelDisp  = true,
      BestLikeStart = 0.98,
      theta       = 1.0
   }
end

-- soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
-- print('soften_length ', calculateEps2(totalBodies, soften_length ))

function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity
    if(run_null_potential == true) then
        print("placing dwarf at origin")
        finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
    else 
        finalPosition, finalVelocity = reverseOrbit{
            potential = potential,
            position  = lbrToCartesian(ctx, Vector.create(l, b, r)),
            velocity  = Vector.create(vx, vy, vz),
            tstop     = revOrbTime,
            dt        = ctx.timestep / 10.0
            }
    end
    
    if(print_reverse_orbit == true) then
        local placeholderPos, placeholderVel = PrintReverseOrbit{
            potential = potential,
            position  = lbrToCartesian(ctx, Vector.create(l, b, r)),
            velocity  = Vector.create(vx, vy, vz),
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
            nbody       = model1Bodies,
            prng        = prng,
            position    = finalPosition,
            velocity    = finalVelocity,
            mass        = dwarfMass,
            scaleRadius = r0,
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
     lambdaStart = -150,
     lambdaEnd = 150,
     lambdaBins = 50,
     
     betaStart = -15,
     betaEnd = 15,
     betaBins = 1
}
end


arg = { ... } -- -- TAKING USER INPUT
assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed") -- STILL EXPECTING SEED AS INPUT FOR THE FUTURE
argSeed = 34086709 -- -- SETTING SEED TO FIXED VALUE
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

if(print_out_parameters) then
    print('forward time=', evolveTime, '\nrev time=',  revOrbTime)
    print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
    print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
    print('total mass solar= ', (mass_d + mass_l) * 222288.47)
    print('rl = ', rscale_l, 'rd = ', rscale_d)
end