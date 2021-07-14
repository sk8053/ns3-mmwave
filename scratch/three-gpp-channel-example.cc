/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2019 SIGNET Lab, Department of Information Engineering,
 * University of Padova
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation;
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/**
* This example shows how to configure the 3GPP channel model classes to
* compute the SNR between two nodes.
* The simulation involves two static nodes which are placed at a certain
* distance from each other and communicates through a wireless channel at
* 2 GHz with a bandwidth of 18 MHz. The default propagation environment is
* 3D-urban macro (UMa) and it can be configured changing the value of the
* string "scenario".
* Each node hosts a SimpleNetDevice and has an antenna array with 4 elements.
*/

#include "ns3/core-module.h"
#include "ns3/three-gpp-channel-model.h"
#include "ns3/three-gpp-antenna-array-model.h"
#include <fstream>
#include "ns3/three-gpp-spectrum-propagation-loss-model.h"
#include "ns3/net-device.h"
#include "ns3/simple-net-device.h"
#include "ns3/node.h"
#include "ns3/node-container.h"
#include "ns3/mobility-model.h"
#include "ns3/constant-position-mobility-model.h"
#include "ns3/lte-spectrum-value-helper.h"
#include "ns3/mmwave-spectrum-value-helper.h"
#include "ns3/channel-condition-model.h"
#include "ns3/three-gpp-propagation-loss-model.h"
#include "ns3/mmwave-helper.h"
#include "fstream"
#include "iostream"
#include <string>
NS_LOG_COMPONENT_DEFINE ("ThreeGppChannelExample");

using namespace ns3;

static Ptr<ThreeGppPropagationLossModel> m_propagationLossModel; //!< the PropagationLossModel object
static Ptr<ThreeGppSpectrumPropagationLossModel> m_spectrumLossModel; //!< the SpectrumPropagationLossModel object

/**
 * Perform the beamforming using the DFT beamforming method
 * \param thisDevice the device performing the beamforming
 * \param thisAntenna the antenna object associated to thisDevice
 * \param otherDevice the device towards which point the beam
 */
static void
DoBeamforming (Ptr<NetDevice> thisDevice, Ptr<ThreeGppAntennaArrayModel> thisAntenna, Ptr<NetDevice> otherDevice)
{
  ThreeGppAntennaArrayModel::ComplexVector antennaWeights;

  // retrieve the position of the two devices
  Vector aPos = thisDevice->GetNode ()->GetObject<MobilityModel> ()->GetPosition ();
  Vector bPos = otherDevice->GetNode ()->GetObject<MobilityModel> ()->GetPosition ();

  // compute the azimuth and the elevation angles
  Angles completeAngle (bPos,aPos);

  double hAngleRadian = fmod (completeAngle.phi, 2.0 * M_PI); // the azimuth angle
  if (hAngleRadian < 0)
  {
    hAngleRadian += 2.0 * M_PI;     
  } 
  double vAngleRadian = completeAngle.theta; // the elevation angle

  // retrieve the number of antenna elements
  int totNoArrayElements = thisAntenna->GetNumberOfElements ();

  // the total power is divided equally among the antenna elements
  double power = 1 / sqrt (totNoArrayElements);

  // compute the antenna weights
  for (int ind = 0; ind < totNoArrayElements; ind++)
    {
      Vector loc = thisAntenna->GetElementLocation (ind);
      double phase = -2 * M_PI * (sin (vAngleRadian) * cos (hAngleRadian) * loc.x
                                  + sin (vAngleRadian) * sin (hAngleRadian) * loc.y
                                  + cos (vAngleRadian) * loc.z);
      antennaWeights.push_back (exp (std::complex<double> (0, phase)) * power);
    }

  // store the antenna weights
  thisAntenna->SetBeamformingVector (antennaWeights);
}

/**
 * Compute the average SNR
 * \param txMob the tx mobility model
 * \param rxMob the rx mobility model
 * \param txPow the transmitting power in dBm
 * \param noiseFigure the noise figure in dB
 */
static void
ComputeSnr (Ptr<MobilityModel> txMob, Ptr<MobilityModel> rxMob, double txPow, double noiseFigure)
{
  // Create the tx PSD using the LteSpectrumValueHelper
  // 100 RBs corresponds to 18 MHz (1 RB = 180 kHz)
  // EARFCN 100 corresponds to 2125.00 MHz
  std::vector<int> activeRbs0 (100);
  for (int i = 0; i < 100 ; i++)
  {
    activeRbs0[i] = i;
  }
  Ptr<SpectrumValue> txPsd = LteSpectrumValueHelper::CreateTxPowerSpectralDensity (2100, 100, txPow, activeRbs0);
  Ptr<SpectrumValue> rxPsd = txPsd->Copy ();
  NS_LOG_DEBUG ("Average tx power " << 10*log10(Sum (*txPsd) * 180e3) << " dB");

  // create the noise PSD
  Ptr<SpectrumValue> noisePsd = LteSpectrumValueHelper::CreateNoisePowerSpectralDensity (2100, 100, noiseFigure);
  NS_LOG_DEBUG ("Average noise power " << 10*log10 (Sum (*noisePsd) * 180e3) << " dB");

  // apply the pathloss
  double propagationGainDb = m_propagationLossModel->CalcRxPower (0, txMob, rxMob);
  NS_LOG_DEBUG ("Pathloss " << -propagationGainDb << " dB");
  double propagationGainLinear = std::pow (10.0, (propagationGainDb) / 10.0);
  *(rxPsd) *= propagationGainLinear;

  NS_LOG_DEBUG (10*log10((*(rxPsd))[0]*180e3));
  // apply the fast fading and the beamforming gain
  rxPsd = m_spectrumLossModel->CalcRxPowerSpectralDensity (rxPsd, txMob, rxMob);
  NS_LOG_DEBUG ("Average rx power " << 10*log10 (Sum (*rxPsd) * 180e3) << " dB");

  // compute the SNR
  NS_LOG_DEBUG ("Average SNR " << 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd)) << " dB");

  // print the SNR and pathloss values in the snr-trace.txt file
  std::ofstream f;
  f.open ("snr-trace.txt", std::ios::out | std::ios::app);
  f << Simulator::Now ().GetSeconds () << " " << 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd)) << " " << propagationGainDb << std::endl;
  f.close ();
}

static void write_file_cluster_only(double *data, std::ofstream *f,  int n_cluster, std::string name){
  //std::ofstream f;
  //f.open ("ray_tracing.txt",  std::ios::app);
  *f << name<<'\t';
  for (uint8_t i = 0 ; i < n_cluster; i++)
        *f<< *(data+i)<<'\t';  
  
    *f << std::endl;


}

/*static void write_file(double *data, std::ofstream *f,  bool angle, int n_cluster, int n_rays, std::string name){
  //std::ofstream f;
  //f.open ("ray_tracing.txt",  std::ios::app);
  uint8_t k = 0;
  *f << name<<'\t';
  for (uint8_t i = 0 ; i < n_cluster; i++){
    
    for (uint8_t j = 0; j<n_rays;j++){
      if (angle == true)
        *f<< *(data+k)<<'\t';
      else
        *f<< *(data+i)<<'\t';
      k++;
    }
    }
    *f << std::endl;
    *f << std::endl;

}*/



static void RayTracing ( Ptr<MobilityModel> * MobPair, double* delay, double * power, double* zod, double*zoa, double* aod, double* aoa, int n_cluster, int n_rays)
{
  double  Pathloss_gain[n_cluster];
  double Zod_c[n_cluster], Zoa_c[n_cluster], Aod_c[n_cluster], Aoa_c[n_cluster];

  //double Zod[n_cluster][n_rays], Zoa[n_cluster][n_rays], Aod[n_cluster][n_rays], Aoa[n_cluster][n_rays];

  Ptr<MobilityModel> txMob  = *MobPair;
  Ptr<MobilityModel> rxMob = *(MobPair+1);
  double propagationGainDb = -1* m_propagationLossModel->CalcRxPower (0, txMob, rxMob);
    Ptr<ChannelConditionModel> cond_model = m_propagationLossModel->GetChannelConditionModel();

  // for loops for priting data and computing power by adding propgation loss
  //uint8_t k = 0;
  for (uint8_t i = 0 ; i < n_cluster; i++){
 
    Pathloss_gain[i] = -10*log10(*(power+i)) + propagationGainDb;
   // for (uint8_t j = 0; j<n_rays;j++){
   //   Zod[i][j] = *(zod+k);
   //   Zoa[i][j] = *(zoa+k);
    //  Aod[i][j] = *(aod+k);
    //  Aoa[i][j] = *(aoa+k);
    //  k++;
     // }
      // sampling one ray per cluster
      Zod_c[i] = *(zod+i);
      Zoa_c[i] = *(zoa+i);
      Aod_c[i] = *(aoa+i);
      Aoa_c[i] = *(aod+i);

  }

  std::ofstream f;
  f.open ("ray_tracing.txt",  std::ios::app);
  f<<"path_number" <<'\t'<< n_cluster  << std::endl;
 uint16_t link_state = cond_model->GetChannelCondition(txMob,rxMob)->GetLosCondition();
 if (link_state == 1)
	link_state = 2; 
  f<<"link_state" <<'\t'<< link_state << std::endl;
  //std::cout<<cond_model->GetChannelCondition(txMob,rxMob)->GetLosCondition();
  f<< "TX"<<'\t'<<txMob->GetPosition().x << '\t' << txMob->GetPosition().y << '\t'<< txMob->GetPosition().z<< std::endl;

  f<< "RX" <<'\t'<<rxMob->GetPosition().x << '\t' << rxMob->GetPosition().y << '\t'<< rxMob->GetPosition().z<< std::endl;
  
  /*write_file(delay, &f, false, n_cluster, n_rays, "delay");
  write_file(Pathloss_gain, &f, false, n_cluster, n_rays, "pathloss");
  write_file(zod, &f, true, n_cluster, n_rays, "zod");
  write_file(zoa, &f, true, n_cluster, n_rays, "zoa");
  write_file(aod, &f, true, n_cluster, n_rays, "aod");
  write_file(aoa, &f, true, n_cluster, n_rays, "aoa");*/

   write_file_cluster_only(delay, &f, n_cluster,  "delay");
  write_file_cluster_only(Pathloss_gain, &f,  n_cluster, "pathloss");
  write_file_cluster_only(Zod_c, &f, n_cluster,  "zod");
  write_file_cluster_only(Zoa_c, &f, n_cluster,  "zoa");
  write_file_cluster_only(Aod_c, &f,  n_cluster, "aod");
  write_file_cluster_only(Aoa_c, &f, n_cluster,  "aoa");

//  NS_LOG_INFO ("sample data"<<"\t" <<Simulator::Now().GetSeconds()<<'\t'<< Delay[0] <<'\t'<< Pathloss_gain[0]<<"\t"<<Aoa[0][0]<<"\t"<<Aod[0][0]<<"\t"<<Zoa[0][0]<<"\t"<<Zod[0][0]);
}

void SetPosition(Ptr<MobilityModel> txMob, Ptr<MobilityModel> rxMob, std::vector<Vector> tx_loc, std::vector<Vector> rx_loc, uint16_t i){
  
  txMob->SetPosition (Vector (tx_loc[i].x,tx_loc[i].y,tx_loc[i].z));
  rxMob->SetPosition (Vector (rx_loc[i].x,rx_loc[i].y,rx_loc[i].z));
 // std::cout << i << std::endl;
  Simulator::Schedule (MilliSeconds (1000), &SetPosition, txMob, rxMob, tx_loc, rx_loc, i+1);
}

int
main (int argc, char *argv[])
{
  std::string line;

  std::vector<Vector> tx_loc, rx_loc;

  std::ifstream f;
  f.open("location.txt");

  while (std::getline(f, line)){
     
      std::istringstream iss(line);
      std::string ind, x, y,z;
      std::getline(iss, ind, '\t' );
      std::getline(iss, x, '\t' );
      std::getline(iss, y, '\t' );
      std::getline(iss, z, '\t' );
    
      if (ind == "tx")
        tx_loc.push_back(Vector(stod(x), stod(y), stod(z)));
      else
        rx_loc.push_back(Vector(stod(x), stod(y), stod(z)));

  }
  f.close();
  //NS_LOG_UNCOND("sizeeee " <<tx_loc.size()) ;
  double frequency = 28.0e9; // operating frequency in Hz (corresponds to EARFCN 2100)
  double txPow = 23.0; // tx power in dBm
  double noiseFigure = 6.0; // noise figure in dB
 // double distance = 10.0; // distance between tx and rx nodes in meters
  uint32_t simTime = 1000*tx_loc.size(); // simulation time in milliseconds
  std::cout <<"simulation time (milliseconds): " << simTime << std::endl; 
  uint32_t timeRes = 10; // time resolution in milliseconds
  std::string scenario = "UMa"; // 3GPP propagation scenario

  Config::SetDefault ("ns3::ThreeGppChannelModel::UpdatePeriod", TimeValue(MilliSeconds (1000))); // update the channel at each iteration
  Config::SetDefault ("ns3::ThreeGppChannelConditionModel::UpdatePeriod", TimeValue(MilliSeconds (10.0))); // do not update the channel condition

  RngSeedManager::SetSeed(1);
  RngSeedManager::SetRun(1);

  // create and configure the factories for the channel condition and propagation loss models
  ObjectFactory propagationLossModelFactory;
  ObjectFactory channelConditionModelFactory;
  if (scenario == "RMa")
  {
    propagationLossModelFactory.SetTypeId (ThreeGppRmaPropagationLossModel::GetTypeId ());
    channelConditionModelFactory.SetTypeId (ThreeGppRmaChannelConditionModel::GetTypeId ());
  }
  else if (scenario == "UMa")
  {
    propagationLossModelFactory.SetTypeId (ThreeGppUmaPropagationLossModel::GetTypeId ());
    channelConditionModelFactory.SetTypeId (ThreeGppUmaChannelConditionModel::GetTypeId ());
  }
  else if (scenario == "UMi-StreetCanyon")
  {
    propagationLossModelFactory.SetTypeId (ThreeGppUmiStreetCanyonPropagationLossModel::GetTypeId ());
    channelConditionModelFactory.SetTypeId (ThreeGppUmiStreetCanyonChannelConditionModel::GetTypeId ());
  }
  else if (scenario == "InH-OfficeOpen")
  {
    propagationLossModelFactory.SetTypeId (ThreeGppIndoorOfficePropagationLossModel::GetTypeId ());
    channelConditionModelFactory.SetTypeId (ThreeGppIndoorOpenOfficeChannelConditionModel::GetTypeId ());
  }
  else if (scenario == "InH-OfficeMixed")
  {
    propagationLossModelFactory.SetTypeId (ThreeGppIndoorOfficePropagationLossModel::GetTypeId ());
    channelConditionModelFactory.SetTypeId (ThreeGppIndoorMixedOfficeChannelConditionModel::GetTypeId ());
  }
  else
  {
    NS_FATAL_ERROR ("Unknown scenario");
  }

  // create the propagation loss model
  m_propagationLossModel = propagationLossModelFactory.Create<ThreeGppPropagationLossModel> ();
  m_propagationLossModel->SetAttribute ("Frequency", DoubleValue (frequency));
  m_propagationLossModel->SetAttribute ("ShadowingEnabled", BooleanValue (true));

  // create the spectrum propagation loss model
  m_spectrumLossModel = CreateObject<ThreeGppSpectrumPropagationLossModel> ();
  m_spectrumLossModel->SetChannelModelAttribute ("Frequency", DoubleValue (frequency));
  m_spectrumLossModel->SetChannelModelAttribute ("Scenario", StringValue (scenario));

  // create the channel condition model and associate it with the spectrum and
  // propagation loss model
  Ptr<ChannelConditionModel> condModel = channelConditionModelFactory.Create<ThreeGppChannelConditionModel> ();
  m_spectrumLossModel->SetChannelModelAttribute ("ChannelConditionModel", PointerValue (condModel));
  m_propagationLossModel->SetChannelConditionModel (condModel);
  //Ptr<MmWaveHelper> mmwaveHelper = CreateObject<MmWaveHelper> ();

  // create the tx and rx nodes
  NodeContainer nodes;
  nodes.Create (2);

  // create the tx and rx devices
 // NetDeviceContainer enbMmWaveDevs = mmwaveHelper->InstallEnbDevice (nodes.Get(0));
 // NetDeviceContainer ueMmWaveDevs = mmwaveHelper->InstallUeDevice (nodes.Get(1));
  //enbNetDev = enbMmWaveDevs.Get (0);
  //ueNetDev = ueMmWaveDevs.Get (0);
  

  Ptr<SimpleNetDevice> txDev = CreateObject<SimpleNetDevice> ();
  Ptr<SimpleNetDevice> rxDev = CreateObject<SimpleNetDevice> ();

  // associate the nodes and the devices
  nodes.Get (0)->AddDevice (txDev);
  txDev->SetNode (nodes.Get (0));
  nodes.Get (1)->AddDevice (rxDev);
  rxDev->SetNode (nodes.Get (1));

  // create the tx and rx mobility models, set the positions
  Ptr<MobilityModel> txMob = CreateObject<ConstantPositionMobilityModel> ();
  txMob->SetPosition (Vector (20.0,0.0,20.0));
  Ptr<MobilityModel> rxMob = CreateObject<ConstantPositionMobilityModel> ();
  rxMob->SetPosition (Vector (0,0.0,1.6));


  // assign the mobility models to the nodes
  nodes.Get (0)->AggregateObject (txMob);
  nodes.Get (1)->AggregateObject (rxMob);

  // create the antenna objects and set their dimensions
  Ptr<ThreeGppAntennaArrayModel> txAntenna = CreateObjectWithAttributes<ThreeGppAntennaArrayModel> ("NumColumns", UintegerValue (2), "NumRows", UintegerValue (2));
  Ptr<ThreeGppAntennaArrayModel> rxAntenna = CreateObjectWithAttributes<ThreeGppAntennaArrayModel> ("NumColumns", UintegerValue (2), "NumRows", UintegerValue (2));

  // Ptr<MmWaveEnbNetDevice> enbNetDevice = StaticCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0));
   //Ptr<ThreeGppAntennaArrayModel> enbAntenna = enbNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
   //enbAntenna->SetAttribute ("IsotropicElements", BooleanValue (true));

  //Ptr<MmWaveEnbNetDevice> enbNetDevice = StaticCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0));
  //Ptr<ThreeGppAntennaArrayModel> ueAntenna = ueNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
   //ueAntenna->SetAttribute ("IsotropicElements", BooleanValue (true));

  // initialize the devices in the ThreeGppSpectrumPropagationLossModel
  m_spectrumLossModel->AddDevice (txDev, txAntenna);
  m_spectrumLossModel->AddDevice (rxDev, rxAntenna);
  auto m_channel =  m_spectrumLossModel->GetChannelModel();
  //std::cout << m_channel->GetChannel(txMob, rxMob,txAntenna, rxAntenna); 
  Ptr<MobilityModel> Mob_pair[2];
  Mob_pair [0] = txMob;
  Mob_pair [1] = rxMob;
   m_channel->TraceConnectWithoutContext ("RayTracing",MakeBoundCallback (&RayTracing, Mob_pair));
  // set the beamforming vectors
  DoBeamforming (txDev, txAntenna, rxDev);
  DoBeamforming (rxDev, rxAntenna, txDev);
  Simulator::Schedule (MilliSeconds (0), &SetPosition,  txMob, rxMob, tx_loc, rx_loc,0);
  
  for (int i = 0; i < floor (simTime / timeRes); i++)
  {
    
    Simulator::Schedule (MilliSeconds (timeRes*i), &ComputeSnr, txMob, rxMob, txPow, noiseFigure);
  }
  Simulator::Stop (MilliSeconds (simTime));
  Simulator::Run ();
  Simulator::Destroy ();
  return 0;
}
