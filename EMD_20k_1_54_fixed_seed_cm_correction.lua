arg = { ... }

assert(#arg == 6, "Expected 6 arguments")
assert(argSeed ~= nil, "Expected seed")
-- argSeed = 34086709
prng = DSFMT.create(argSeed)
print("seed: ", argSeed)
-- npa = new parameter arrangement --

evolveTime       = tonumber(arg[1])
reverseOrbitTime = tonumber(arg[1]) / tonumber(arg[2])
rscale_l         = tonumber(arg[3])
light_r_ratio    = tonumber(arg[4])
mass_l           = tonumber(arg[5])
light_mass_ratio = tonumber(arg[6])

model1Bodies = 20
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
print("evolve time: ", evolveTime)
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

-- This function takes the table of bodies produced and zeroes the center of mass 
-- and center of momentum. It then shifts the center of mass and center of momentum
-- to the expected value for its position in the orbit.
function cm_correction(firstModel, dwarf_starting_position, dwarf_starting_velocity)

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
        print('old positions: ' , firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z)
        print('old velocities: ', firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
    end
    cm_x  = cm_x  / dwarfMass
    cm_y  = cm_y  / dwarfMass
    cm_z  = cm_z  / dwarfMass
    cm_vx = cm_vx / dwarfMass
    cm_vy = cm_vy / dwarfMass
    cm_vz = cm_vz / dwarfMass
    print('old CM: ', cm_x, cm_y, cm_z)
    print('old CMV: ', cm_vx, cm_vy, cm_vz)
    for i, v in ipairs(firstModel)
    do
        x_new = firstModel[i].position.x - cm_x + dwarf_starting_position.x
        y_new = firstModel[i].position.y - cm_y + dwarf_starting_position.y
        z_new = firstModel[i].position.z - cm_z + dwarf_starting_position.z
        
        vx_new = firstModel[i].velocity.x - cm_vx + dwarf_starting_velocity.x
        vy_new = firstModel[i].velocity.y - cm_vy + dwarf_starting_velocity.y
        vz_new = firstModel[i].velocity.z - cm_vz + dwarf_starting_velocity.z
        
        firstModel[i].position = Vector.create(x_new, y_new, z_new)
        firstModel[i].velocity = Vector.create(vx_new, vy_new, vz_new)
    end
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
        print('new positions: ' , firstModel[i].position.x, firstModel[i].position.y, firstModel[i].position.z)
        print('new velocities: ', firstModel[i].velocity.x, firstModel[i].velocity.y, firstModel[i].velocity.z)
    end
    cm_x  = cm_x  / dwarfMass
    cm_y  = cm_y  / dwarfMass
    cm_z  = cm_z  / dwarfMass
    cm_vx = cm_vx / dwarfMass
    cm_vy = cm_vy / dwarfMass
    cm_vz = cm_vz / dwarfMass
    print('new CM: ', cm_x, cm_y, cm_z)
    print('new CMV: ', cm_vx, cm_vy, cm_vz)
  return firstModel
end

-- Also required
-- position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6))
-- velocity  = Vector.create(-156, 79, 107),
function makeBodies(ctx, potential)
  local firstModel
  --this puts a particle at the current postion and velocity of the Orphan stream
  --calculates the reverse orbit along the potential and finds the required
  --center of mass position and velocity.
  local dwarf_starting_position, dwarf_starting_velocity = reverseOrbit{
      potential = potential,
      position  = lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)),
      velocity  = Vector.create(-156, 79, 107),
      tstop     = reverseOrbitTime,
      dt        = ctx.timestep / 10.0,
  }
--     print('final position: ', dwarf_starting_position.x, dwarf_starting_position.y, dwarf_starting_position.z)
--     print('final velocity: ',dwarf_starting_velocity.x, dwarf_starting_velocity.y, dwarf_starting_velocity.z)
--     print("position: ", lbrToCartesian(ctx, Vector.create(218, 53.5, 28.6)))
--     print("velocity: ", Vector.create(-156, 79, 107))
  firstModel = predefinedModels.isotropic{
      nbody       = model1Bodies,
      prng        = prng,
      position    = dwarf_starting_position,
      velocity    = dwarf_starting_velocity,
      mass1       = mass_l,
      mass2       = mass_d,
      scaleRadius1 = rscale_l,
      scaleRadius2 = rscale_d,
      ignore      = true
  }
  
--   cm_correction(firstModel, dwarf_starting_position, dwarf_starting_velocity)

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
     lambdaBins = 50,
     betaStart = -40,
     betaEnd = 40,
     betaBins = 1
}
end