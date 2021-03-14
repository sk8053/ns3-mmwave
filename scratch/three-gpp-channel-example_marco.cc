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
#include "ns3/mmwave-spectrum-value-helper.h"
#include "ns3/lte-spectrum-value-helper.h"
#include <ns3/mmwave-phy-mac-common.h>
#include <ns3/mmwave-helper.h>
#include "ns3/channel-condition-model.h"
#include "ns3/three-gpp-propagation-loss-model.h"

NS_LOG_COMPONENT_DEFINE ("ThreeGppChannelExample");

using namespace ns3;
using namespace mmwave;

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
  //NS_LOG_UNCOND("Position A is " << aPos);
  //NS_LOG_UNCOND("Position B is " << bPos);

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
  //NS_LOG_UNCOND("Tot number of elements is " << totNoArrayElements);

  // the total power is divided equally among the antenna elements
  // double power = 1 / sqrt (totNoArrayElements);
  //NS_LOG_UNCOND("Power is " << 10*log10(power));

  // compute the antenna weights
  for (int ind = 0; ind < totNoArrayElements; ind++)
    {
      Vector loc = thisAntenna->GetElementLocation (ind);
      //NS_LOG_UNCOND("Location of the antenna: " << loc);
      double phase = -2 * M_PI * (sin (vAngleRadian) * cos (hAngleRadian) * loc.x
                                  + sin (vAngleRadian) * sin (hAngleRadian) * loc.y
                                  + cos (vAngleRadian) * loc.z);
      // MM: remove power
      antennaWeights.push_back (exp (std::complex<double> (0, phase)));
      //antennaWeights.push_back (exp (std::complex<double> (0, phase)) * power);
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
ComputeSnr (Ptr<MobilityModel> txMob, Ptr<MobilityModel> rxMob, Ptr<MmWaveEnbNetDevice> enbDev, Ptr<MmWavePhyMacCommon> phyMac, double txPow, double noiseFigure)
{
	int n_rbs = phyMac->GetNumChunks();
	std::vector<int> activeRbs0 (n_rbs);
	for (int i = 0; i < n_rbs ; i++)
	{
		activeRbs0[i] = i;
	}

  Ptr<SpectrumValue> txPsd = MmWaveSpectrumValueHelper::CreateTxPowerSpectralDensity (phyMac, txPow, activeRbs0);
  Ptr<SpectrumValue> rxPsd = txPsd->Copy ();
  //NS_LOG_UNCOND("Total TX power across all chunks: " << 10*log10(Sum (*txPsd) * phyMac->GetChunkWidth()) << " dB");
  //NS_LOG_UNCOND("TX power on chunk 1: " << 10*log10((*txPsd)[0] * phyMac->GetChunkWidth()) << " dB");

  // create the noise PSD
  Ptr<SpectrumValue> noisePsd = MmWaveSpectrumValueHelper::CreateNoisePowerSpectralDensity (phyMac, noiseFigure);
  //NS_LOG_UNCOND("Noise power across all chunks: " << 10*log10 (Sum (*noisePsd) * phyMac->GetChunkWidth()) << " dB");
  //NS_LOG_UNCOND("Noise power on chunk 1: " << 10*log10((*noisePsd)[0] * phyMac->GetChunkWidth()) << " dB");

  // apply the pathloss
  double propagationGainDb = m_propagationLossModel->CalcRxPower (0, txMob, rxMob);
  NS_LOG_UNCOND("Pathloss " << -propagationGainDb << " dB");
  double propagationGainLinear = std::pow (10.0, (propagationGainDb) / 10.0);
  *(rxPsd) *= propagationGainLinear;

  // apply the fast fading and the beamforming gain
  rxPsd = m_spectrumLossModel->CalcRxPowerSpectralDensity (rxPsd, txMob, rxMob);
  NS_LOG_UNCOND("Average rx power " << 10*log10 (Sum (*rxPsd) * phyMac->GetChunkWidth()) << " dB");

  // compute the SNR
  NS_LOG_UNCOND("Average SNR " << 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd)) << " dB");

  // print the SNR and pathloss values in the snr-trace.txt file
  std::ofstream f;
  f.open ("snr-trace.txt", std::ios::out | std::ios::app);
  f << Simulator::Now ().GetSeconds () << " " << 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd)) << " " << propagationGainDb << std::endl;
  f.close ();
}

static void RayTracing ( Ptr<MobilityModel> * MobPair, double* delay, double * power, double* zod, double*zoa, double* aod, double* aoa, int n_cluster, int n_rays)

{
  double Delay[n_cluster], Pathloss_gain[n_cluster];
  //double Zod_c[n_cluster], Zoa_c[n_cluster], Aod_c[n_cluster], Aoa_c[n_cluster];
  double Zod[n_cluster][n_rays], Zoa[n_cluster][n_rays], Aod[n_cluster][n_rays], Aoa[n_cluster][n_rays];

  Ptr<MobilityModel> txMob  = *MobPair;
  Ptr<MobilityModel> rxMob = *(MobPair+1);

  double propagationGainDb = -1* m_propagationLossModel->CalcRxPower (0, txMob, rxMob);

  // for loops for printing data and computing power by adding the propagation loss
  uint8_t k = 0;
  for (uint8_t i = 0 ; i < n_cluster; i++){
    Delay[i] = *(delay+i);
    Pathloss_gain[i] = -10*log10(*(power+i)) + propagationGainDb;
    for (uint8_t j = 0; j<n_rays;j++){
      Zod[i][j] = *(zod+k);
      Zoa[i][j] = *(zoa+k);
      Aod[i][j] = *(aod+k);
      Aoa[i][j] = *(aoa+k);
      k++;
      }
      // sampling one ray per cluster
    //  Zod_c[i] = Zod[i][0];
    //  Zoa_c[i] = Zoa[i][0];
    //  Aod_c[i] = Aod[i][0];
    //  Aoa_c[i] = Aoa[i][0];
  }

  NS_LOG_UNCOND ("sample data"<<"\t" <<Simulator::Now().GetSeconds()<<'\t'<< Delay[0] <<'\t'<< Pathloss_gain[0]<<"\t"<<Aoa[0][0]<<"\t"<<Aod[0][0]<<"\t"<<Zoa[0][0]<<"\t"<<Zod[0][0]);
}

int
main (int argc, char *argv[])
{

  double txPow = 35.0; // tx power in dBm
  double noiseFigure = 9.0; // noise figure in dB
  double distance = 100.0; // distance between tx and rx nodes in meters
  uint32_t simTime = 10000; // simulation time in milliseconds
  uint32_t timeRes = 100; // time resolution in milliseconds
  uint16_t enbAntennaNum = 8; // The number of antenna elements on one axis
  uint16_t ueAntennaNum = 4; // The number of antenna elements on one axis
  std::string scenario = "UMi-StreetCanyon"; // 3GPP propagation scenario

  Config::SetDefault ("ns3::ThreeGppChannelModel::UpdatePeriod", TimeValue(MilliSeconds (1))); // update the channel at each iteration
  Config::SetDefault ("ns3::ThreeGppChannelConditionModel::UpdatePeriod", TimeValue(MilliSeconds (0.0))); // do not update the channel condition

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

  	// create the tx and rx nodes
    NodeContainer nodes;
    nodes.Create (2);

    // create the tx and rx mobility models, set the positions
    Ptr<MobilityModel> txMob = CreateObject<ConstantPositionMobilityModel> ();
    txMob->SetPosition (Vector (0.0,0.0,5.0));
    Ptr<MobilityModel> rxMob = CreateObject<ConstantPositionMobilityModel> ();
    rxMob->SetPosition (Vector (distance,0.0,1.6));

    // assign the mobility models to the nodes
    nodes.Get (0)->AggregateObject (txMob);
    nodes.Get (1)->AggregateObject (rxMob);


    // Create the MmWave helper
    Ptr<MmWaveHelper> mmwaveHelper = CreateObject<MmWaveHelper> ();

    // Create the tx and rx devices
    NetDeviceContainer enbMmWaveDevs = mmwaveHelper->InstallEnbDevice (nodes.Get (0));
    NetDeviceContainer ueMmWaveDevs = mmwaveHelper->InstallUeDevice (nodes.Get (1));
    Ptr<MmWaveEnbNetDevice> enbNetDevice = StaticCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0));
    Ptr<MmWaveUeNetDevice> ueNetDevice = StaticCast<MmWaveUeNetDevice> (ueMmWaveDevs.Get (0));

    // associate the nodes and the devices
	nodes.Get (0)->AddDevice (enbNetDevice);
	enbNetDevice->SetNode (nodes.Get (0));
	nodes.Get (1)->AddDevice (ueNetDevice);
	ueNetDevice->SetNode (nodes.Get (1));

	Ptr<ThreeGppAntennaArrayModel> enbAntenna = enbNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
	enbAntenna->SetAttribute ("NumRows", UintegerValue (enbAntennaNum));
	enbAntenna->SetAttribute ("NumColumns", UintegerValue (enbAntennaNum));
	enbAntenna->SetAttribute ("DowntiltAngle", DoubleValue (DegreesToRadians (12.0)));
	enbAntenna->SetAttribute ("BearingAngle" , DoubleValue (DegreesToRadians (0.0)));
	enbAntenna->SetAttribute ("ElementGain" , DoubleValue (8.0));

	Ptr<ThreeGppAntennaArrayModel> ueAntenna = ueNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
	ueAntenna->SetAttribute ("NumRows", UintegerValue (ueAntennaNum));
	ueAntenna->SetAttribute ("NumColumns", UintegerValue (ueAntennaNum));
	//ueAntenna->SetAttribute ("DowntiltAngle" , DoubleValue (DegreesToRadians (180.0)));
	ueAntenna->SetAttribute ("BearingAngle" , DoubleValue (DegreesToRadians (0.0)));
	ueAntenna->SetAttribute ("ElementGain" , DoubleValue (8.0));

	// Num3: 277 RBs corresponds to 400 MHz (1 RB = 1.44 MHz - SCS 120 KHz, 12 subcarriers)
    Ptr<MmWavePhyMacCommon> phyMac = DynamicCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0))->GetPhy ()->GetConfigurationParameters ();
    //enbNetDevice->GetPhy ()->GetConfigurationParameters ();
	//phyMac->SetNumerology (MmWavePhyMacCommon::Numerology::NrNumerology3);

	// create the 3GPP propagation loss model
	m_propagationLossModel = propagationLossModelFactory.Create<ThreeGppPropagationLossModel> ();
	m_propagationLossModel->SetAttribute ("Frequency", DoubleValue (phyMac->GetCenterFrequency()));
	m_propagationLossModel->SetAttribute ("ShadowingEnabled", BooleanValue (false));

	// create the 3GPP spectrum propagation loss model
	m_spectrumLossModel = CreateObject<ThreeGppSpectrumPropagationLossModel> ();
	m_spectrumLossModel->SetChannelModelAttribute ("Frequency", DoubleValue (phyMac->GetCenterFrequency()));
	m_spectrumLossModel->SetChannelModelAttribute ("Scenario", StringValue (scenario));

	// create the 3GPP channel condition model
	Ptr<ChannelConditionModel> condModel = channelConditionModelFactory.Create<ThreeGppChannelConditionModel> ();
	m_spectrumLossModel->SetChannelModelAttribute ("ChannelConditionModel", PointerValue (condModel));
	m_propagationLossModel->SetChannelConditionModel (condModel);

	// initialize the devices in the ThreeGppSpectrumPropagationLossModel
	m_spectrumLossModel->AddDevice (enbNetDevice, enbAntenna);
	m_spectrumLossModel->AddDevice (ueNetDevice, ueAntenna);
 // mmwaveHelper->AttachToClosestEnb(ueNetDevice, enbNetDevice);
  
    auto m_channel =  m_spectrumLossModel->GetChannelModel();

    Ptr<MobilityModel> Mob_pair[2];
    Mob_pair [0] = txMob;
    Mob_pair [1] = rxMob;

    m_channel->TraceConnectWithoutContext ("RayTracing", MakeBoundCallback (&RayTracing, Mob_pair));

    // set the beamforming vectors
    DoBeamforming (enbNetDevice, enbAntenna, ueNetDevice);
    DoBeamforming (ueNetDevice, ueAntenna, enbNetDevice);

  for (int i = 0; i < floor (simTime / timeRes); i++)
  {
	  Simulator::Schedule (MilliSeconds (timeRes*i), &ComputeSnr, txMob, rxMob, enbNetDevice, phyMac, txPow, noiseFigure);
  }

  Simulator::Run ();
  Simulator::Destroy ();
  return 0;
}
