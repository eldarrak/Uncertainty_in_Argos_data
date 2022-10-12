##############
# propagate argos errors
#############
# ver 0.1
# 12_10_22
# ER

#############################
# PART 1 - uncertainty
#############################
library(dplyr)
library(move)
library(ctmm)
#setwd('./Argos_ellipsoids')

cur_track_move<-readRDS('App-Output_ Workflow_Instance_001__Movebank__2022-10-11_08-15-33.rds') %>% subset(argos.lc %in% c("1","2","3")) %>%  subset(!is.na(location.long)) %>% subset(argos.semi.major!=0) %>% subset(argos.semi.major!=0) %>% subset(location.lat>52.8) 
cur_track_ctmm<-as.telemetry(cur_track_move)

# here we can also, download not a movebank file, but anything else..
# but I suggest we start from the move object..

# automated guestimate for calibrated data
GUESS <- ctmm.guess(cur_track_ctmm, CTMM=ctmm(error=TRUE),interactive=FALSE)
FIT <- ctmm.select(cur_track_ctmm,GUESS,trace=TRUE,cores=6)

n.simulations=100  
 
orig.df<-as.data.frame(as.data.frame(cur_track_ctmm), precompute=TRUE)
# I will make a brick with dimnesions (coord, time, simID)

sim.brick<-array(NA, c( nrow(orig.df), 2, n.simulations))

Sim_1<-simulate(FIT, data=cur_track_ctmm, precompute=TRUE, t=cur_track_ctmm$t) 

   for (i in 1:n.simulations) {
      cat('\r', i)
      Sim_1<-simulate(FIT, data=cur_track_ctmm, precompute=-1, t=cur_track_ctmm$t) # now we could sample 
      sim.brick[,,i]<- suppressWarnings(SpatialPoints(coords = cbind(Sim_1@.Data[[2]], Sim_1@.Data[[3]]), proj4string = CRS(FIT@info$projection))) %>% 
      spTransform(crs('+proj=longlat +datum=WGS84 +no_defs'))    %>% coordinates()
	  }

# Now we just want to save and then overlay with the map, when needed.

# Now I want to overlay this output with the RUth_map..
Output.sim<-cur_track_move
attr(Output.sim,"sim.error") <- sim.brick

saveRDS(Output.sim, file='track.with.100.simulated.errors.RDS')

#########################################
# Part 2. Spatial overlay example
#########################################
# ok, now we want to overlay with Ruth map...
library(rgeos)

## Load Habitat map
Cur_sp<-readOGR('./Ruth_map/Std_Dev_polygons_rdnew.shp')

#Make a buffer around the map (will need it to distinguish between the points falling beyond the map range and those falling on 'not assigned' habitat within the map)
Buffer_50<-gBuffer(Cur_sp,  width=50) 

all_out<-c()

for (i in 1:dim(attr(Output.sim,"sim.error"))[3]) {
      cat('\r', i)
      Sim_SPts<-suppressWarnings(SpatialPoints(coords = attr(Output.sim,"sim.error")[,,i], proj4string = CRS(Output.sim@proj4string@projargs))) 
				Sim_SPts<-spTransform(Sim_SPts,crs(Cur_sp))
				#remove points outside the area
				overlay_area<-extract(Cur_sp, Sim_SPts)
				overlay_buffer<-extract(Buffer_50,Sim_SPts)
				names(overlay_buffer)<-c("point.ID.buffer","poly.ID.buffer" )
				habitat_NA_within<-cbind(overlay_area,overlay_buffer)
				habitat_NA_within$habitat[is.na(habitat_NA_within$habitat) & !is.na(habitat_NA_within$poly.ID.buffer)]<-'not_assigned'

	  all_out<-cbind(all_out, habitat_NA_within$habitat)
   }

   Habitat_types<-as.data.frame(t(apply(all_out, 1, FUN=function(x) unlist(table(factor(x,  levels=c('high_change', 'intermediate', 'low_change', 'not_assigned' )), useNA='always')))))
   names(Habitat_types)<-paste0('ctmm_simulation_', as.character(names(Habitat_types)))

   cur_track_move$most.present.habitat_ctmm_sim<-
   apply(Habitat_types, 1, FUN=function(x) c('high_change', 'intermediate', 'low_change','not_assigned')[which.max(x[1:4])])
   cur_track_move<-cbind(cur_track_move, Habitat_types)  