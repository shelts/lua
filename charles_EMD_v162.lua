arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
prng = DSFMT.create(argSeed)

evolveTime       = tonumber(arg[1])
revOrbTime       = tonumber(arg[2])
rscale_l         = tonumber(arg[3])
rscale_d         = tonumber(arg[4])
mass_l           = tonumber(arg[5])
mass_d           = tonumber(arg[6])

mass_l = mass_l / 222288.47
mass_d = mass_d / 222288.47
model1Bodies = 20000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.62"
run_null_potential = false
print_reverse_orbit = false

print('forward time=', evolveTime, '\nrev time=',  revOrbTime)
print('mass_l sim=', mass_l, '\nmass_d sim=', mass_d)
print('light mass solar=', mass_l * 222288.47, '\ndark mass solar=', mass_d * 222288.47)
print('total mass solar= ', (mass_d + mass_l) * 222288.47)
print('rl = ', rscale_l, 'rd = ', rscale_d)

function makePotential()
   if(run_null_potential == true) then
       print("running in null potential")
       return nil
   else
        return  Potential.create{
            spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
            disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
            halo      = Halo.nfw{ vhalo = 73, scaleLength = 22.250}
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
    print('timestep ', t, tmp)
    
    return t
end


function makeContext()
   soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = get_timestep(),
      eps2       = calculateEps2(totalBodies, soften_length ),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
print('soften_length ', calculateEps2(totalBodies, soften_length ))

function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity
  if(run_null_potential == true) then
      print("placing dwarf at origin")
      finalPosition, finalVelocity = Vector.create(0, 0, 0), Vector.create(0, 0, 0)
  else 
    finalPosition, finalVelocity = reverseOrbit{
        potential = potential,
        position  = lbrToCartesian(ctx, Vector.create(45, 47.20, 11.51)),
        velocity  = Vector.create(-114.87, 166.06, 57.55),
        tstop     = revOrbTime,
        dt        = ctx.timestep / 10.0
        }
  end
    
  if(print_reverse_orbit == true) then
    local placeholderPos, placeholderVel = PrintReverseOrbit{
        potential = potential,
        position  = lbrToCartesian(ctx, Vector.create(45, 47.20, 11.51)),
        velocity  = Vector.create(-114.87, 166.06, 57.55),
        tstop     = .14,
        tstopf    = .20,
        dt        = ctx.timestep / 10.0
    }
    print('Printing reverse orbit')
  end
  
print(lbrToCartesian(ctx, Vector.create(45, 46.93, 11.87)), Vector.create(-122.78, 157.32, 64.90))
print(finalPosition, finalVelocity)

--   firstModel = predefinedModels.isotropic{
--       nbody       = model1Bodies,
--       prng        = prng,
--       position    = finalPosition,
--       velocity    = finalVelocity,
--       mass1       = mass_l,
--       mass2       = mass_d,
--       scaleRadius1 = rscale_l,
--       scaleRadius2 = rscale_d,
--       ignore      = true
--   }
  
  
     firstModel = predefinedModels.plummer{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass        = mass_d,
      scaleRadius = rscale_d,
      ignore      = false
   }
  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -150,
     lambdaEnd = 150,
     lambdaBins = 100,
     betaStart = -100,
     betaEnd = 100,
     betaBins = 1
}
end


