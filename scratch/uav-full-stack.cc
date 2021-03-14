/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2020 SIGNET Lab, Department of Information Engineering,
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
 * This example shows how to configure a full stack simulation using the
 * QdChannelModel and traces from Remcom.
 * The simulation involves the following scenario:
 * UE --------------------------> UE
 *               BS 
 */

#include <fstream>
#include "ns3/core-module.h"
#include "ns3/three-gpp-spectrum-propagation-loss-model.h"
#include "ns3/simple-net-device.h"
#include "ns3/node-container.h"
#include "ns3/constant-position-mobility-model.h"
#include "ns3/qd-channel-model.h"
#include "ns3/mmwave-helper.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/config-store.h"
#include "ns3/mmwave-point-to-point-epc-helper.h"

NS_LOG_COMPONENT_DEFINE ("RemcomModelExample");

using namespace ns3;
using namespace mmwave;

// main variables
Ptr<QdChannelModel> qdModel;
Ptr<PacketSink> sinkApp;              // Sink Application at the UE
Ptr<OutputStreamWrapper> stream1, stream2, stream3; // Output stream 
Ptr<SpectrumValue> txPsd, noisePsd;
Ptr<NetDevice> enbNetDev, ueNetDev;
uint32_t timeRes = 10;                // Time resolution in ms
uint32_t lastRxBytes = 0;

static void 
Rx (Ptr<const Packet> packet, const Address &from)
{
  NS_LOG_DEBUG (Simulator::Now ().GetSeconds () << "\t" << packet->GetSize ());
  //*stream1->GetStream () << Simulator::Now ().GetSeconds () << "\t" << packet->GetSize () << std::endl;
}

static void
ComputeE2eThroughput ()
{
  uint32_t totRxBytes = sinkApp->GetTotalRx ();   // Total rx bytes 
  uint32_t rxBytes = totRxBytes - lastRxBytes;
  lastRxBytes = totRxBytes;
  double thr = rxBytes * 8.0 / (timeRes * 1e-3) / 1e6;  // Throughput in Mbps
  NS_LOG_UNCOND (Simulator::Now ().GetSeconds () << "\t" << thr << " Mbps");
  *stream2->GetStream () << Simulator::Now ().GetSeconds () << "\t" << thr << std::endl;

  Simulator::Schedule (MilliSeconds (timeRes), &ComputeE2eThroughput);
}

static void
Sinr (uint64_t imsi, SpectrumValue& oldSinr, SpectrumValue& newSinr)
{
  double sinr = Sum (newSinr) / (newSinr.GetSpectrumModel ()->GetNumBands ());
  NS_LOG_DEBUG (Simulator::Now ().GetSeconds () << "\t" << 10*log10 (sinr) << " dB");
  //*stream3->GetStream () << Simulator::Now ().GetSeconds () << "\t" << sinr << std::endl;
}

static void
ComputeSnr ()
{
  Ptr<MmWaveSpectrumPhy> enbSpectrumPhy = DynamicCast<MmWaveEnbNetDevice> (enbNetDev)->GetPhy ()->GetDlSpectrumPhy ();
  Ptr<MmWaveSpectrumPhy> ueSpectrumPhy = DynamicCast<MmWaveUeNetDevice> (ueNetDev)->GetPhy ()->GetDlSpectrumPhy ();
  enbSpectrumPhy->ConfigureBeamforming (ueNetDev);

  /*Ptr<ThreeGppAntennaArrayModel> enbAntenna = enbNetDev->GetNode ()->GetObject<ThreeGppAntennaArrayModel> ();
  NS_LOG_UNCOND (enbAntenna);
  Ptr<ThreeGppAntennaArrayModel> ueAntenna = ueNetDev->GetNode ()->GetObject<ThreeGppAntennaArrayModel> ();
  auto enbVector = enbAntenna->GetBeamformingVector ();
  auto ueVector = ueAntenna->GetBeamformingVector ();
  for (uint32_t i = 0; i < enbVector.size (); ++i)
  {
    enbVector.at (i) = std::polar (1.0,0.0);
  }
  for (uint32_t i = 0; i < ueVector.size (); ++i)
  {
    ueVector.at (i) = std::polar (1.0,0.0);
  }
  enbAntenna->SetBeamformingVector (enbVector);
  ueAntenna->SetBeamformingVector (ueVector);*/

  // retrieve the PropagationLossModel and the SpectrumPropagationLossModel
  Ptr<SpectrumChannel> sc = enbSpectrumPhy->GetSpectrumChannel ();
  Ptr<SpectrumPropagationLossModel> splm = sc->GetSpectrumPropagationLossModel ();
  
  // compute the sinr
  Ptr<MobilityModel> enbMob = enbNetDev->GetNode ()->GetObject<MobilityModel> ();
  Ptr<MobilityModel> ueMob = ueNetDev->GetNode ()->GetObject<MobilityModel> ();
  
  Ptr<SpectrumValue> rxPsd = splm->CalcRxPowerSpectralDensity (txPsd, enbMob, ueMob);
  //NS_LOG_UNCOND (rxPsd->GetValuesN ());
  double snr = 10 * log10 (Sum (*rxPsd) / Sum (*noisePsd));
  NS_LOG_UNCOND ("Snr [dB]= " << snr);
  *stream3->GetStream () << Simulator::Now ().GetSeconds () << "\t" << Sum (*rxPsd) / Sum (*noisePsd) << std::endl;

}

int
main (int argc, char *argv[])
{
  std::string qdFilesPath = "contrib/qd-channel/model/QD/"; // The path of the folder with the QD scenarios
  std::string scenario = "Uav"; // The name of the scenario
  uint32_t interPacketInterval = 1e3; // App inter packet arrival [us]
  double txPower = 23.0; // Transmitted power for both eNB and UE [dBm]
  double noiseFigure = 6.0; // Noise figure for both eNB and UE [dB]
  uint16_t enbAntennaNum = 1; // The number of antenna elements for the gNBs antenna arrays, assuming a square architecture
  uint16_t ueAntennaNum = 1; // The number of antenna elements for the UE antenna arrays, assuming a square architecture
  uint32_t appPacketSize = 1460; // Application packet size [B]
  uint32_t bandwidth = 400e6;    // Bandwidth of the system
  bool fullStack = false;

  CommandLine cmd;
  cmd.AddValue ("qdFilesPath", "The path of the folder with the QD scenarios", qdFilesPath);
  cmd.AddValue ("scenario", "The name of the scenario", scenario);
  cmd.AddValue ("ipi", "App inter packet arrival [us]", interPacketInterval);
  cmd.AddValue ("txPower", "Transmitted power for both eNB and UE [dBm]", txPower);
  cmd.AddValue ("noiseFigure", "Noise figure for both eNB and UE [dB]", noiseFigure);
  cmd.AddValue ("enbAntennaNum", "The number of antenna elements for the gNBs antenna arrays, assuming a square architecture", enbAntennaNum);
  cmd.AddValue ("ueAntennaNum", "The number of antenna elements for the UE antenna arrays, assuming a square architecture", ueAntennaNum);
  cmd.AddValue ("appPacketSize", "Application packet size [B]", appPacketSize);
  cmd.Parse (argc, argv);

  // Setup
  LogComponentEnableAll (LOG_PREFIX_ALL);

  bool harqEnabled = true;
  bool rlcAmEnabled = true;

  Config::SetDefault ("ns3::MmWavePhyMacCommon::Bandwidth", DoubleValue (bandwidth));
  Config::SetDefault ("ns3::MmWaveEnbPhy::TxPower", DoubleValue (txPower));
  Config::SetDefault ("ns3::MmWaveEnbPhy::NoiseFigure", DoubleValue (noiseFigure));
  Config::SetDefault ("ns3::MmWaveUePhy::TxPower", DoubleValue (txPower));
  Config::SetDefault ("ns3::MmWaveUePhy::NoiseFigure", DoubleValue (noiseFigure));

  Config::SetDefault ("ns3::MmWaveHelper::RlcAmEnabled", BooleanValue (rlcAmEnabled));
  Config::SetDefault ("ns3::MmWaveHelper::HarqEnabled", BooleanValue (harqEnabled));
  Config::SetDefault ("ns3::MmWaveFlexTtiMacScheduler::HarqEnabled", BooleanValue (harqEnabled));

  // Create the tx and rx nodes
  NodeContainer ueNodes;
  NodeContainer enbNodes;
  enbNodes.Create (1);
  ueNodes.Create (1);

  // Create the tx and rx mobility models, set the positions to be equal to the
  // initial positions of the nodes in the ray tracer
  Ptr<MobilityModel> ueRefMob = CreateObject<ConstantPositionMobilityModel> ();
  ueRefMob->SetPosition (Vector (0, 0, 30));
  Ptr<MobilityModel> enb1Mob = CreateObject<ConstantPositionMobilityModel> ();
  enb1Mob->SetPosition (Vector (0, 0, 10));

  // Assign the mobility models to the nodes
  enbNodes.Get (0)->AggregateObject (enb1Mob);
  ueNodes.Get (0)->AggregateObject (ueRefMob);

  // Configure the channel
  Config::SetDefault ("ns3::MmWaveHelper::PathlossModel", StringValue (""));
  Config::SetDefault ("ns3::MmWaveHelper::ChannelModel", StringValue ("ns3::ThreeGppSpectrumPropagationLossModel"));
  qdModel = CreateObject<QdChannelModel> (qdFilesPath, scenario);
  Time simTime = qdModel->GetQdSimTime ();
  Config::SetDefault ("ns3::ThreeGppSpectrumPropagationLossModel::ChannelModel", PointerValue (qdModel));

  // Create the MmWave helper
  Ptr<MmWaveHelper> mmwaveHelper = CreateObject<MmWaveHelper> ();
  mmwaveHelper->SetMmWaveEnbNetDeviceAttribute("AntennaNum", UintegerValue (enbAntennaNum));
  mmwaveHelper->SetMmWaveUeNetDeviceAttribute("AntennaNum", UintegerValue (ueAntennaNum));

  mmwaveHelper->SetSchedulerType ("ns3::MmWaveFlexTtiMacScheduler");
  Ptr<MmWavePointToPointEpcHelper>  epcHelper = CreateObject<MmWavePointToPointEpcHelper> ();
  mmwaveHelper->SetEpcHelper (epcHelper);
  mmwaveHelper->SetHarqEnabled (harqEnabled);

  // Create a single RemoteHost
  Ptr<Node> pgw = epcHelper->GetPgwNode ();
  NodeContainer remoteHostContainer;
  remoteHostContainer.Create (1);
  Ptr<Node> remoteHost = remoteHostContainer.Get (0);
  InternetStackHelper internet;
  internet.Install (remoteHostContainer);

  // Create the Internet
  PointToPointHelper p2ph;
  p2ph.SetDeviceAttribute ("DataRate", DataRateValue (DataRate ("100Gb/s")));
  p2ph.SetDeviceAttribute ("Mtu", UintegerValue (1500));
  p2ph.SetChannelAttribute ("Delay", TimeValue (Seconds (0.010)));
  NetDeviceContainer internetDevices = p2ph.Install (pgw, remoteHost);
  Ipv4AddressHelper ipv4h;
  ipv4h.SetBase ("1.0.0.0", "255.0.0.0");
  Ipv4InterfaceContainer internetIpIfaces = ipv4h.Assign (internetDevices);
  // Interface 0 is localhost, 1 is the p2p device
  Ipv4StaticRoutingHelper ipv4RoutingHelper;
  Ptr<Ipv4StaticRouting> remoteHostStaticRouting = ipv4RoutingHelper.GetStaticRouting (remoteHost->GetObject<Ipv4> ());
  remoteHostStaticRouting->AddNetworkRouteTo (Ipv4Address ("7.0.0.0"), Ipv4Mask ("255.0.0.0"), 1);

  // Create the tx and rx devices
  NetDeviceContainer enbMmWaveDevs = mmwaveHelper->InstallEnbDevice (enbNodes);
  NetDeviceContainer ueMmWaveDevs = mmwaveHelper->InstallUeDevice (ueNodes);
  ueNetDev = ueMmWaveDevs.Get (0);
  enbNetDev = enbMmWaveDevs.Get (0);
  AsciiTraceHelper asciiTraceHelper;

  Ptr<MmWaveEnbNetDevice> enbNetDevice = StaticCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0));
  Ptr<ThreeGppAntennaArrayModel> enbAntenna = enbNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
  //enbAntenna->SetAttribute ("DowntiltAngle", DoubleValue (DegreesToRadians (12.0)));
  // enbAntenna->SetAttribute ("BearingAngle" , DoubleValue (DegreesToRadians (0.0)));
  //enbAntenna->SetAttribute ("ElementGain" , DoubleValue (8.0));
  enbAntenna->SetAttribute ("IsotropicElements", BooleanValue (true));
  
  Ptr<MmWaveUeNetDevice> ueNetDevice = StaticCast<MmWaveUeNetDevice> (ueMmWaveDevs.Get (0));
  Ptr<ThreeGppAntennaArrayModel> ueAntenna = ueNetDevice->GetPhy ()->GetDlSpectrumPhy ()->GetBeamformingModel ()->GetAntenna ();
  //ueAntenna->SetAttribute ("DowntiltAngle" , DoubleValue (DegreesToRadians (180.0)));
  //ueAntenna->SetAttribute ("ElementGain" , DoubleValue (8.0));
  ueAntenna->SetAttribute ("IsotropicElements", BooleanValue (true));

  if (!fullStack)
    {
      // create the psd of the transmitted signal and of the noise
      stream3 = asciiTraceHelper.CreateFileStream ("snr-vs-time.txt");
      Ptr<MmWavePhyMacCommon> mwpmc = DynamicCast<MmWaveEnbNetDevice> (enbMmWaveDevs.Get (0))->GetPhy ()->GetConfigurationParameters ();
      std::vector <int> activeRbs;
      for (uint32_t i = 0; i < mwpmc->GetNumChunks (); i++)
      {
        activeRbs.push_back (i);
      }
      txPsd = MmWaveSpectrumValueHelper::CreateTxPowerSpectralDensity (mwpmc, txPower, activeRbs);
      noisePsd = MmWaveSpectrumValueHelper::CreateNoisePowerSpectralDensity (mwpmc, noiseFigure);
      
      Time timeStep = MilliSeconds (timeRes);
      for (auto i = 0; i < simTime.GetMicroSeconds () / timeStep.GetMicroSeconds (); i++)
      {
        Simulator::Schedule (timeStep * i, &ComputeSnr);
      }
    }
  else
    {
      // Install the IP stack on the UEs
      internet.Install (ueNodes);
      Ipv4InterfaceContainer ueIpIface;
      ueIpIface = epcHelper->AssignUeIpv4Address (NetDeviceContainer (ueMmWaveDevs));
      // Assign IP address to UEs, and install applications
      Ptr<Node> ueNode = ueNodes.Get (0);
      // Set the default gateway for the UE
      Ptr<Ipv4StaticRouting> ueStaticRouting = ipv4RoutingHelper.GetStaticRouting (ueNode->GetObject<Ipv4> ());
      ueStaticRouting->SetDefaultRoute (epcHelper->GetUeDefaultGatewayAddress (), 1);

      // This performs the attachment of each UE to a specific eNB
      mmwaveHelper->AttachToEnbWithIndex (ueMmWaveDevs.Get (0), enbMmWaveDevs, 0);

      // Add apps
      uint16_t dlPort = 1234;
      uint16_t ulPort = 2000;
      uint16_t otherPort = 3000;
      ApplicationContainer clientApps;
      ApplicationContainer serverApps;
      ++ulPort;
      ++otherPort;
      PacketSinkHelper dlPacketSinkHelper ("ns3::UdpSocketFactory", InetSocketAddress (Ipv4Address::GetAny (), dlPort));
      serverApps.Add (dlPacketSinkHelper.Install (remoteHost));

      /*UdpClientHelper dlClient (ueIpIface.GetAddress (0), dlPort);
      dlClient.SetAttribute ("Interval", TimeValue (MicroSeconds (interPacketInterval)));
      dlClient.SetAttribute ("MaxPackets", UintegerValue (1000000));
      dlClient.SetAttribute ("PacketSize", UintegerValue (appPacketSize));
      clientApps.Add (dlClient.Install (remoteHost));*/

      OnOffHelper onOffClient ("ns3::UdpSocketFactory", InetSocketAddress (internetIpIfaces.GetAddress (1), dlPort));
      onOffClient.SetAttribute ("PacketSize", UintegerValue (appPacketSize));
      onOffClient.SetAttribute ("DataRate", DataRateValue (DataRate ("1000Mbps")));
      onOffClient.SetAttribute ("OnTime", StringValue ("ns3::ConstantRandomVariable[Constant=1.0e6]"));
      onOffClient.SetAttribute ("OffTime", StringValue ("ns3::ConstantRandomVariable[Constant=0]"));
      clientApps.Add (onOffClient.Install (ueNodes.Get (0)));

      serverApps.Start (Seconds (0.001));
      clientApps.Start (Seconds (0.001));
      mmwaveHelper->EnableTraces ();

      stream1 = asciiTraceHelper.CreateFileStream ("rx-packet-trace.txt");
      stream2 = asciiTraceHelper.CreateFileStream ("thr-vs-time.txt");
      sinkApp = StaticCast<PacketSink> (serverApps.Get (0));
      sinkApp->TraceConnectWithoutContext ("Rx", MakeCallback (&Rx));

      Ptr<MmWaveUePhy> uePhy = StaticCast <MmWaveUeNetDevice> (ueMmWaveDevs.Get (0))->GetPhy ();
      uePhy->TraceConnectWithoutContext ("ReportCurrentCellRsrpSinr", MakeCallback (&Sinr));

      Simulator::Schedule (MilliSeconds (timeRes), &ComputeE2eThroughput);
    }

  Simulator::Stop (simTime);
  Simulator::Run ();
  Simulator::Destroy ();
  return 0;
}
