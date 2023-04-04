#include "TopTrackerConstruction.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4RegionStore.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"

#include "Identifier/TtID.h"

#include <TString.h>

#include <iomanip>
#include <exception>
using namespace CLHEP;
DECLARE_TOOL(TopTrackerConstruction);

TopTrackerConstruction::TopTrackerConstruction(const std::string& name)
    : ToolBase(name)
{
    wall_mapping = NULL;
    logicAirTT = 0;
    initVariables();
   
}

TopTrackerConstruction::~TopTrackerConstruction() {
  if(wall_mapping != NULL) delete wall_mapping;
  
}

G4LogicalVolume* 
TopTrackerConstruction::getLV() {
  // if already initialized
  if (logicAirTT) {
    return logicAirTT;
  }
  initMaterials();
  
  makeAirTTLogical();
  
  makeWallLogical();
  makePlaneLogical();
  makePanelLogical();
  makePanelTapeLogical();
  makeCoatingLogical();
  makeStripLogical();
  
  makeWallPhysical();
  for(std::set<int>::iterator iPiWMN = panel_in_wall_mask_needed.begin(); iPiWMN != panel_in_wall_mask_needed.end(); ++iPiWMN){
    makePlanePhysical(*iPiWMN);
  }
  for(std::set<int>::iterator iPiPMN = panel_in_plane_mask_needed.begin(); iPiPMN != panel_in_plane_mask_needed.end(); ++iPiPMN){
    makePanelPhysical(*iPiPMN);
  }
  makePanelTapePhysical();
  makeCoatingPhysical();
  makeStripPhysical();
  
  return logicAirTT;
}

bool
TopTrackerConstruction::inject(std::string motherName, IDetElement* other, IDetElementPos* /*pos*/) {
  // Get the mother volume in current DetElem.
  G4LogicalVolume* mothervol = 0;
  if ( motherName == "lAirTT" ) {
    mothervol = logicAirTT;
  }
  if (not mothervol) {
    // don't find the volume.
    return false;
  }
  
  // retrieve the daughter's LV
  G4LogicalVolume* daughtervol = other->getLV();
  
  if (not daughtervol) {
    return false;
  }
  
  
  return true;
}

double 
TopTrackerConstruction::geom_info(const std::string& param, const int wallnumber)
{
  if (param == "WallNum") {
    return std::next(m_wall_posx.begin(),wallnumber)->first;
  } else if (param == "Wall.X") {
    return m_wall_posx[wallnumber];
  } else if (param == "Wall.Y") {
    return m_wall_posy[wallnumber];
  } else if (param == "Wall.Z") {
    return m_wall_posz[wallnumber];
  } else if (param == "Wall.Layer") {
    return m_wall_layer[wallnumber];
  }  else if (param == "Wall.Column") {
    return m_wall_column[wallnumber];
  }  else if (param == "Wall.Row") {
    return m_wall_row[wallnumber];
  } else if (param == "Panel.XY") {
    return m_panel_posxy[wallnumber];
  } else if (param == "Bar.XY") {
    return m_bar_posxy[wallnumber];
  } else if (param == "Plane.Z") {
    return m_plane_posz[wallnumber];
  }
  
  else {
    // don't recognize, throw exception
    throw std::invalid_argument("Unknown param: "+param);
  }
}

void TopTrackerConstruction::readWallPlacement(std::string & txt_file){
  if(wall_mapping != NULL) delete wall_mapping;
  wall_mapping = new std::map<int, struct wall_pos>;
  std::ifstream infile(txt_file);
  while(true){
    int wallID;
    struct wall_pos this_wall;
    std::string line;
    if(!std::getline(infile, line)) break;
    if(infile.eof()) break;
    if(line[0] == '#' || line[0] == '\n'){
      // comment or empty line!
      continue;
    }

    std::istringstream line_ss(line);
    line_ss >> wallID
      >> this_wall.xc >> this_wall.yc >> this_wall.zc
      >> this_wall.angle
      >> std::setbase(0) >> this_wall.panel_mask ;
    if(line_ss.fail()) continue;

    (*wall_mapping)[wallID] = this_wall;

    panel_in_wall_mask_needed .insert(this_wall.panel_mask);
    panel_in_plane_mask_needed.insert(this_wall.panel_mask & 0xF);
    panel_in_plane_mask_needed.insert((this_wall.panel_mask>>4) & 0xF);

  }
  infile.close();
}

void
TopTrackerConstruction::initVariables() {
  
  
  m_box_x= 48./2.*m;
  m_box_y= 48./2.*m;
  m_box_z= 8.40/2.*m;
  
  m_chimney_R=50*cm;
  m_chimney_Z=2.0*m;
  // the designed chimney radius is 40cm, the chimney height above the water is ~3.5m.
  
  m_lengthBar= 686/2.*cm;
  m_thicknessBar= 1/2.*cm;
  m_widthBar= 2.6/2.*cm;
  m_spaceBar=0.01*cm;
  
  m_lengthCoating= 686/2.*cm;
  m_thicknessCoating= 1.03/2.*cm;
  m_widthCoating= 2.63/2.*cm;
  
  m_lengthPanelTape=686/2.*cm;
  m_thicknessPanelTape=1.21/2.*cm;
  m_widthPanelTape=169.13/2.*cm;
  
  m_lengthPanel=686.12/2.*cm;
  m_thicknessPanel=1.33/2.*cm;
  m_widthPanel=169.25/2.*cm;
  m_spacePanel=0.01*cm;
  
  m_lengthPlane=686.12/2.*cm;
  m_thicknessPlane=1.33/2.*cm;
  m_widthPlane=677.03/2.*cm;
  m_vspacePlane=0.1*cm;
  
  m_wall_x=686.12/2.*cm;
  m_wall_y=686.12/2.*cm;
  m_wall_z= 2.76/2.*cm;
  
  m_floor_level = -m_box_z;

  std::string wall_placement_file  = std::getenv("TOPTRACKERROOT");
  wall_placement_file += "/include/TopTrackerWallPlacement.txt";
  readWallPlacement(wall_placement_file);
  
  //IsOverlap = true;
  IsOverlap = false;
  
}

void 
TopTrackerConstruction::initMaterials() {
  Aluminium = G4Material::GetMaterial("Aluminium");
  Scintillator = G4Material::GetMaterial("Scintillator");
  Adhesive = G4Material::GetMaterial("Adhesive");
  TiO2Coating = G4Material::GetMaterial("TiO2Coating");
  Air = G4Material::GetMaterial("Air");
  
}

void
TopTrackerConstruction::makeAirTTLogical() {
  BoxsolidAirTT = new G4Box("BoxsAirTT",
			    m_box_x,
			    m_box_y,
			    m_box_z);
  
  solidChimney = new G4Tubs("Cylinder",
			    0,
			    m_chimney_R,
			    m_chimney_Z,
			    0,
			    twopi); 
  
  G4ThreeVector zTrans(0, 0,-m_box_z + m_chimney_Z);
  
  solidAirTT =  new G4SubtractionSolid("sAirTT",
				       BoxsolidAirTT,
				       solidChimney,
				       0,
				       zTrans);
  
  
  logicAirTT = new G4LogicalVolume(solidAirTT,
				   Air,
				   "lAirTT",
				   0,
				   0,
				   0);
}

void
TopTrackerConstruction::makeStripLogical() {
  
  solidBar = new G4Box("sBar",
		       m_lengthBar,
		       m_widthBar,
		       m_thicknessBar);
  
  logicBar = new G4LogicalVolume(solidBar,
				 Scintillator,
				 "lBar",
				 0,
				 0,
				 0);
  
}
void
TopTrackerConstruction::makeStripPhysical() {
  
  
  physiBar = new G4PVPlacement(0,              // no rotation
			       G4ThreeVector(0,0,0), // at (x,y,z)
			       logicBar,    // its logical volume
			       "pBar",       // its name
			       logicCoating,  // its mother  volume
			       false,           // no boolean operations
			       0,
			       IsOverlap);              // no particular field
  
}
void
TopTrackerConstruction::makeCoatingLogical() {
  
  solidCoating = new G4Box("sBar",
			   m_lengthCoating,
			   m_widthCoating,
			   m_thicknessCoating);
  
  logicCoating = new G4LogicalVolume(solidCoating,
				     TiO2Coating,
				     "lCoating",
				     0,
				     0,
				     0);
  
}
void
TopTrackerConstruction::makeCoatingPhysical() {
  
  char coat_name[30];
  double y0(0);
  double yy(0);
  
  y0=-(64*2*m_widthCoating + 63*m_spaceBar)/2. + m_widthCoating;
  
  for(int jj=0;jj<64;jj++) {
    
    sprintf(coat_name,"pCoating_%02d_",jj);
    
    yy=y0+jj * (m_spaceBar + 2*m_widthCoating);
    
    physiCoating[jj] = new G4PVPlacement(0,              // no rotation
					 G4ThreeVector(0,yy,0), // at (x,y,z)
					 logicCoating,    // its logical volume
					 coat_name,       // its name
					 logicPanelTape,  // its mother  volume
					 false,           // no boolean operations
					 jj,  //copy number
					 IsOverlap);
    
    m_bar_posxy[jj]=yy;
    
    
  }
  
}

void
TopTrackerConstruction::makePanelTapeLogical() {
  
  solidPanelTape = new G4Box("sPanelTape",
			      m_lengthPanelTape,
			      m_widthPanelTape,
			      m_thicknessPanelTape);
  
  logicPanelTape = new G4LogicalVolume(solidPanelTape,
					Adhesive,
					"lPanelTape",
					0,
					0,
					0);
  
}
void
TopTrackerConstruction::makePanelTapePhysical() {
  
  
  physiPanelTape = new G4PVPlacement(0,              // no rotation
				      G4ThreeVector(0,0,0), // at (x,y,z)
				      logicPanelTape,    // its logical volume
				      "pPanelTape",       // its name
				      logicPanel,  // its mother  volume
				      false,           // no boolean operations
				      0,
				      IsOverlap);              // no particular field
  
}

void
TopTrackerConstruction::makePanelLogical() {
  
  solidPanel = new G4Box("sPanel",
			  m_lengthPanel,
			  m_widthPanel,
			  m_thicknessPanel);
  
  logicPanel = new G4LogicalVolume(solidPanel,
				    Aluminium,
				    "lPanel",
				    0,
				    0,
				    0);
  
}

void
TopTrackerConstruction::makePanelPhysical(int panel_mask) {
  
  static int copyno = 0;
  char panel_name[30];
  double y0(0);
  double yy(0);
  
  y0=-(4*2*m_widthPanel + 3*m_spacePanel)/2. + m_widthPanel;
  
  for(int jj=0;jj<4;jj++) {
    if(((panel_mask>>jj) & 0x1) == 0) continue;
    
    sprintf(panel_name,"pPanel_%d_%x_",jj,panel_mask);
    
    yy=y0+jj * (m_spacePanel + 2*m_widthPanel);
    
    physiPanel.push_back(new G4PVPlacement(0,              // no rotation
					G4ThreeVector(0,yy,0), // at (x,y,z)
					logicPanel,    // its logical volume
					panel_name,       // its name
					logicPlane[panel_mask],  // its mother  volume
					false,           // no boolean operations
					copyno++,  //copy number
					IsOverlap));
       
    m_panel_posxy[jj]=yy;
    
  }
  
}
void
TopTrackerConstruction::makePlaneLogical() {
  
  solidPlane = new G4Box("sPlane",
                         m_lengthPlane,
                         m_widthPlane,
                         m_thicknessPlane);
  
  for(std::set<int>::iterator iPiPMN = panel_in_plane_mask_needed.begin(); iPiPMN != panel_in_plane_mask_needed.end(); ++iPiPMN){
    TString logic_plane_name = TString::Format("lPlane%x_", *iPiPMN);
    logicPlane[*iPiPMN] = new G4LogicalVolume(solidPlane,
                                              Air,
                                              logic_plane_name.Data(),
                                              0,
                                              0,
                                              0);
  }
  
}

void
TopTrackerConstruction::makePlanePhysical(int plane_mask) {
  static int copyno = 0;
  char plane_name[30];
  double z0(0);
  double zz(0);
  
  z0=-(2*2*m_thicknessPlane + m_vspacePlane)/2. + m_thicknessPlane;
  
  G4RotationMatrix MRot;
  MRot.rotateZ(-M_PI/2.*rad); 
  
  for(int jj=0;jj<2;jj++) {
    int panel_mask = (plane_mask>>(4*jj)) & 0xF ;
    if(panel_mask == 0) continue;
    
    sprintf(plane_name,"pPlane_%d_%02x_",jj,plane_mask);
    
    zz=z0+jj * (m_vspacePlane + 2*m_thicknessPlane);
    if(jj==0) //bottom plane, "vertical walls"
      {
	physiPlane.push_back(new G4PVPlacement(
					   G4Transform3D(MRot,G4ThreeVector(0,0,zz)),
					   logicPlane[panel_mask],    // its logical volume
					   plane_name,       // its name
					   logicWall[plane_mask],  // its mother  volume
					   false,           // no boolean operations
					   copyno++,  //copy number
					   IsOverlap));
      }
    else  //upper plane, "horizontal planes"
      {
	physiPlane.push_back(new G4PVPlacement(
					   0,G4ThreeVector(0,0,zz),
					   logicPlane[panel_mask],    // its logical volume
					   plane_name,       // its name
					   logicWall[plane_mask],  // its mother  volume
					   false,           // no boolean operations
					   copyno++,  //copy number
					   IsOverlap));
      }
    
    m_plane_posz[jj]=zz;
  }
  
}

void
TopTrackerConstruction::makeWallLogical() {
  
  solidWall = new G4Box("sWall",
			m_wall_x,
			m_wall_y,
			m_wall_z);
  
  for(std::set<int>::iterator iPiWMN = panel_in_wall_mask_needed.begin(); iPiWMN != panel_in_wall_mask_needed.end(); ++iPiWMN){
    TString logic_wall_name = TString::Format("lWall%02x_", *iPiWMN);
    logicWall[*iPiWMN] = new G4LogicalVolume(solidWall,
                                            Air,
                                            logic_wall_name.Data(),
                                            0,
                                            0,
                                            0);
  }

}

void TopTrackerConstruction::makeWallPhysical() {
  static int copyno = 0;
  if(wall_mapping == NULL || wall_mapping->size() == 0){
    LogError << "Wall Mapping not defined. Aborting." << std::endl;
    throw std::runtime_error("Bad TT Wall Mapping");
  }

  std::map<int, struct wall_pos>::iterator wall_map_iter;
  for(wall_map_iter = wall_mapping->begin(); wall_map_iter != wall_mapping->end(); ++wall_map_iter){
    int wallID = wall_map_iter->first;
    struct wall_pos TT_wp = wall_map_iter->second;

    TString wall_name = TString::Format("pWall_%03d_", wallID);
    G4ThreeVector pos(TT_wp.xc, TT_wp.yc, m_floor_level + TT_wp.zc);
    G4RotationMatrix rot;
    rot.rotateZ(TT_wp.angle * deg);

    physiWall.push_back( new G4PVPlacement(
        G4Transform3D(rot, pos),
        logicWall[TT_wp.panel_mask],
        wall_name.Data(),
        logicAirTT,
        false,
        copyno,
        IsOverlap));

    m_wall_posx  [wallID]=pos.x();
    m_wall_posy  [wallID]=pos.y();
    m_wall_posz  [wallID]=pos.z();
    m_wall_layer [wallID]=TtID::layer (TtID::id(wallID,0));
    m_wall_column[wallID]=TtID::column(TtID::id(wallID,0));
    m_wall_row   [wallID]=TtID::row   (TtID::id(wallID,0));

  }

}


