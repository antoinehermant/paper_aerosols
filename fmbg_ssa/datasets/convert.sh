cdo selvar,dR_sp*,*sra*,albedo,aclcov historical_simple-plumes_ssa093_BOT_PRP_yearmean.nc dR_sp.nc
cdo divc,0.00000055 -selvar,*aod* historical_simple-plumes_ssa093_BOT_PRP_yearmean.nc aod.nc
cdo merge dR_sp.nc aod.nc temp.nc