arg = { ... }

assert(#arg == 7, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
argSeed = 34086709
print("seed: ", argSeed)
prng = DSFMT.create(argSeed)

-- npa = new parameter arrangement --

evolveTime       = tonumber(arg[1])
revOrbTime       = tonumber(arg[1]) / tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])
version          = tonumber(arg[7])

print('version: ', version)
if version == 0 then
    correcting = 'y'
elseif version == 1 then
    correcting = 'n'
end

model1Bodies = 20000
totalBodies = model1Bodies

nbodyLikelihoodMethod = "EMD"
nbodyMinVersion = "1.54"

--fitting ml directly. mass ratio as originally defined
--radius ratio now defined the way mass ratio is.
--bit of a construct because rscale_t is not really a thing
--works better because this method makes r_d more sensitive to changes in radius_ratio

dwarfMass = mass_l / light_mass_ratio
rscale_t  = rscale_l / light_r_ratio
rscale_d  = rscale_t *  (1.0 - light_r_ratio)
mass_d    = dwarfMass * (1.0 - light_mass_ratio)

--print(evolveTime)
-- print(mass_d, mass_l)
function makePotential()
   return  Potential.create{
      spherical = Spherical.spherical{ mass  = 1.52954402e5, scale = 0.7 },
      disk      = Disk.miyamotoNagai{ mass = 4.45865888e5, scaleLength = 6.5, scaleHeight = 0.26 },
      halo      = Halo.logarithmic{ vhalo = 73, scaleLength = 12.0, flattenZ = 1.0 }
   }
end

function get_timestep()
    --Mass of a single dark matter sphere enclosed within light rscale
    mass_enc_d = mass_d * (rscale_l)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    --Mass of a single light matter sphere enclosed within dark rscale
    mass_enc_l = mass_l * (rscale_d)^3 * ( (sqr(rscale_l)+ sqr(rscale_d) ) )^(-3.0/2.0)

    s1 = cube(rscale_l) / (mass_enc_d + mass_l)
    s2 = cube(rscale_d) / (mass_enc_l + mass_d)
    
    --return the smaller time step
    if(s1 < s2) then
        s = s1
    else
        s = s2
    end
    
    -- I did it this way so there was only one place to change the time step. 
    t = (1/100) * sqrt( pi_4_3 * s)
    return t
end


function makeContext()
   soften_length = (mass_l * rscale_l + mass_d  * rscale_d) / (mass_d + mass_l)
   return NBodyCtx.create{
      timeEvolve = evolveTime,
      timestep   = get_timestep(),
      eps2       = calculateEps2(totalBodies, soften_length),
      criterion  = "NewCriterion",
      useQuad    = true,
      theta      = 1.0
   }
end

function cm_correction(firstModel, finalPosition, finalVelocity)
    cm_x = 0.0
    cm_y = 0.0
    cm_z = 0.0
    cm_vx = 0.0
    cm_vy = 0.0
    cm_vz = 0.0
    for i, v in ipairs(firstModel)
    do
        cm_x  = cm_x  + ( firstModel[i].mass * firstModel[i].position.x) 
        cm_y  = cm_y  + ( firstModel[i].mass * firstModel[i].position.y) 
        cm_z  = cm_z  + ( firstModel[i].mass * firstModel[i].position.z) 
        
        cm_vx = cm_vx + ( firstModel[i].mass * firstModel[i].velocity.x) 
        cm_vy = cm_vy + ( firstModel[i].mass * firstModel[i].velocity.y) 
        cm_vz = cm_vz + ( firstModel[i].mass * firstModel[i].velocity.z) 
        if correcting == 'n' then
--             print(firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z, firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
        end
        
--         print('old positions: ', firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z)
--         print('old velocities: ', firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
--         print(firstModel[i].mass)
    end
    cm_x  = cm_x  / dwarfMass
    cm_y  = cm_y  / dwarfMass
    cm_z  = cm_z  / dwarfMass
    cm_vx = cm_vx / dwarfMass
    cm_vy = cm_vy / dwarfMass
    cm_vz = cm_vz / dwarfMass
--     print('old CM: ', cm_x, cm_y, cm_z)
--     print('old CMV: ', cm_vx, cm_vy, cm_vz)
    
    if correcting == 'y' then
--         print("correcting: ")
        for i, v in ipairs(firstModel)
        do
            x_new = firstModel[i].position.x - cm_x + finalPosition.x
            y_new = firstModel[i].position.y - cm_y + finalPosition.y
            z_new = firstModel[i].position.z - cm_z + finalPosition.z
            
            vx_new = firstModel[i].velocity.x - cm_vx + finalVelocity.x
            vy_new = firstModel[i].velocity.y - cm_vy + finalVelocity.y
            vz_new = firstModel[i].velocity.z - cm_vz + finalVelocity.z
            
            firstModel[i].position = Vector.create(x_new, y_new, z_new)
            firstModel[i].velocity = Vector.create(vx_new, vy_new, vz_new)
        end
    end
    cm_x = 0.0
    cm_y = 0.0
    cm_z = 0.0
    cm_vx = 0.0
    cm_vy = 0.0
    cm_vz = 0.0
    t_mass = 0.0
    for i, v in ipairs(firstModel)
    do
        cm_x  = cm_x  + ( firstModel[i].mass * firstModel[i].position.x) 
        cm_y  = cm_y  + ( firstModel[i].mass * firstModel[i].position.y) 
        cm_z  = cm_z  + ( firstModel[i].mass * firstModel[i].position.z) 
        
        cm_vx = cm_vx + ( firstModel[i].mass * firstModel[i].velocity.x) 
        cm_vy = cm_vy + ( firstModel[i].mass * firstModel[i].velocity.y) 
        cm_vz = cm_vz + ( firstModel[i].mass * firstModel[i].velocity.z) 
        if correcting == 'y' then
--             print(firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z, firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
        end
--         print('new positions: ' , firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z)
--         print('new velocities: ', firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
        t_mass = t_mass + firstModel[i].mass
    end
    
    print(t_mass)
    cm_x  = cm_x  / dwarfMass
    cm_y  = cm_y  / dwarfMass
    cm_z  = cm_z  / dwarfMass
    cm_vx = cm_vx / dwarfMass
    cm_vy = cm_vy / dwarfMass
    cm_vz = cm_vz / dwarfMass
--     print('new CM: ', cm_x, cm_y, cm_z)
--     print('new CMV: ', cm_vx, cm_vy, cm_vz)
  return firstModel
end

-- Also required
-- position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6))
function makeBodies(ctx, potential)
  local firstModel
  local finalPosition, finalVelocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity  = Vector.create(-156, 79, 107),
      tstop     = revOrbTime,
      dt        = ctx.timestep / 10.0
  }
    print('final position: ', finalPosition.x, finalPosition.y, finalPosition.z)
    print('final velocity: ',finalVelocity.x, finalVelocity.y, finalVelocity.z)
  firstModel = predefinedModels.isotropic{
      nbody       = model1Bodies,
      prng        = prng,
      position    = finalPosition,
      velocity    = finalVelocity,
      mass1       = mass_l,
      mass2       = mass_d,
      scaleRadius1 = rscale_l,
      scaleRadius2 = rscale_d,
      ignore      = true
  }
  
  cm_correction(firstModel, finalPosition, finalVelocity)

  return firstModel
end

function makeHistogram()
    return HistogramParams.create{
     --Orphan Stream coordinate transformation angles
     phi = 128.79,
     theta = 54.39,
     psi = 90.70,
     lambdaStart = -75,
     lambdaEnd = 50,
     lambdaBins = 50,
     betaStart = -40,
     betaEnd = 40,
     betaBins = 1
}
end


