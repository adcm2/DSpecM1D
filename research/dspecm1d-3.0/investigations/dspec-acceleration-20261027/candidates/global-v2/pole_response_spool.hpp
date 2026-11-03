#pragma once

#include "global_gsh_tile.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace dspec_global_gsh {
struct NativeSide {
  int element=-1,node=-1;
  double radiusMetres=0.0,radiusModel=0.0,widthMetres=0.0;
  bool fluid=false;
};
struct SpoolAudit { std::int32_t frequencyBin=-1; Complex z{},adjointMixing{}; double basisReconstructionRelative=0.0; };
struct FrequencySideRecords { PoleRecord forward,unitAdjoint;L0Record l0;SpoolAudit audit; };

/** Append-only compact pole records. Index tables live on disk; resident memory
 * is bounded by native-side metadata, per-key completion flags, and frequency
 * audits. Use one Reader per worker because its stream position is mutable. */
class PoleResponseSpoolWriter {
public:
  PoleResponseSpoolWriter(const std::string &path,int lmax,std::size_t bins,
      std::vector<NativeSide> sides,std::uint64_t byteCap)
      :path_(path),lmax_(lmax),bins_(bins),sides_(std::move(sides)),byteCap_(byteCap) {
    if(lmax_<1||lmax_>1536||bins_==0||sides_.empty()||byteCap_==0||
       bins_>std::numeric_limits<std::int32_t>::max()+std::size_t{1}||
       bins_>std::numeric_limits<std::size_t>::max()/sides_.size()/2)
      throw std::invalid_argument("invalid pole response spool dimensions or cap");
    const auto parent=std::filesystem::path(path_).parent_path();
    if(!parent.empty())std::filesystem::create_directories(parent);
    if(std::filesystem::exists(path_))throw std::runtime_error("response spool path already exists");
    for(std::size_t i=0;i<sides_.size();++i) {
      const auto &s=sides_[i];
      if(s.element<0||s.node<0||!(std::isfinite(s.radiusMetres)&&s.radiusMetres>0.0&&
         std::isfinite(s.radiusModel)&&s.radiusModel>0.0&&std::isfinite(s.widthMetres)&&s.widthMetres>0.0)||
         !sideIndex_.emplace(std::pair{s.element,std::bit_cast<std::uint64_t>(s.radiusModel)},i).second)
        throw std::invalid_argument("native-side spool metadata is invalid or duplicated");
    }
    const std::uint64_t labels=5ull*static_cast<std::uint64_t>(lmax_+1);
    const std::uint64_t perPole=sizeof(std::uint64_t)+9*sizeof(double)+
        labels*(2*sizeof(std::int32_t)+6*sizeof(Complex));
    if(perPole>(std::numeric_limits<std::uint64_t>::max()-4*sizeof(Complex))/2)
      throw std::overflow_error("response spool record estimate overflows");
    const std::uint64_t perKey=2*perPole+4*sizeof(Complex);
    if(bins_>std::numeric_limits<std::uint64_t>::max()/sides_.size()/perKey)
      throw std::overflow_error("response spool estimate overflows");
    const std::uint64_t payloadEstimate=static_cast<std::uint64_t>(bins_)*sides_.size()*perKey;
    const std::uint64_t indexBytes=static_cast<std::uint64_t>(bins_)*sides_.size()*4*sizeof(std::uint64_t);
    const std::uint64_t fixedBytes=128+indexBytes+static_cast<std::uint64_t>(sides_.size())*37+
        static_cast<std::uint64_t>(bins_)*(sizeof(std::int32_t)+2*sizeof(Complex)+sizeof(double));
    if(payloadEstimate>std::numeric_limits<std::uint64_t>::max()-fixedBytes||
       payloadEstimate+fixedBytes>byteCap_)
      throw std::runtime_error("worst-case compact response spool exceeds configured byte cap");
    const auto space=std::filesystem::space(parent.empty()?std::filesystem::current_path():parent);
    if(payloadEstimate+fixedBytes>space.available)
      throw std::runtime_error("insufficient free disk space for worst-case compact response spool");
    const std::size_t keys=bins_*sides_.size();
    if(keys>512u*1024u*1024u/3)
      throw std::runtime_error("resident response-spool completion metadata exceeds its 512 MiB cap");
    recordSeen_.assign(keys*2,0);l0Seen_.assign(keys,0);
    audits_.assign(bins_,{});auditSet_.assign(bins_,false);
    stream_.open(path_,std::ios::binary|std::ios::in|std::ios::out|std::ios::trunc);
    if(!stream_)throw std::runtime_error("cannot create pole response spool");
    const std::array<char,8> magic{'D','S','G','S','H','S','P','1'};
    write(magic.data(),magic.size());writeValue(std::uint32_t{1});writeValue(static_cast<std::int32_t>(lmax_));
    writeValue(static_cast<std::uint64_t>(bins_));writeValue(static_cast<std::uint64_t>(sides_.size()));
    for(const auto &s:sides_){writeValue(static_cast<std::int32_t>(s.element));writeValue(static_cast<std::int32_t>(s.node));
      writeValue(s.radiusMetres);writeValue(s.radiusModel);writeValue(s.widthMetres);writeValue(static_cast<std::uint8_t>(s.fluid));}
    recordOffsetsOffset_=static_cast<std::uint64_t>(stream_.tellp());
    for(std::size_t i=0;i<keys*2;++i)writeValue(std::uint64_t{0});
    l0OffsetsOffset_=static_cast<std::uint64_t>(stream_.tellp());
    for(std::size_t i=0;i<keys*2;++i)writeValue(std::uint64_t{0});
    auditOffset_=static_cast<std::uint64_t>(stream_.tellp());
    for(std::size_t i=0;i<bins_;++i){writeValue(std::int32_t{-1});writeValue(Complex{});writeValue(Complex{});writeValue(0.0);}
    dataOffset_=static_cast<std::uint64_t>(stream_.tellp());
  }

  void setAudit(std::size_t bin,SpoolAudit audit) {
    requireOpen();if(bin>=bins_||audit.frequencyBin!=static_cast<std::int64_t>(bin)||!finite(audit.z)||!finite(audit.adjointMixing)||
       !std::isfinite(audit.basisReconstructionRelative)||audit.basisReconstructionRelative<0.0)
      throw std::invalid_argument("invalid response-spool adjoint audit");
    if(auditSet_[bin])throw std::runtime_error("duplicate response-spool frequency audit");
    audits_[bin]=audit;auditSet_[bin]=true;
  }

  void append(std::size_t bin,const PoleRecord &record) {
    requireOpen();if(bin>=bins_||record.element<0||!std::isfinite(record.radiusModel))
      throw std::invalid_argument("invalid response-spool pole record");
    const auto sideIt=sideIndex_.find({record.element,std::bit_cast<std::uint64_t>(record.radiusModel)});
    if(sideIt==sideIndex_.end())throw std::runtime_error("pole record does not match a native side");
    const std::size_t side=sideIt->second,key=bin*sides_.size()+side;
    if(record.synthesisSet==std::numeric_limits<std::uint64_t>::max()||
       record.synthesisSet==std::numeric_limits<std::uint64_t>::max()-1) {
      const std::size_t family=record.synthesisSet==std::numeric_limits<std::uint64_t>::max()?0:1;
      if(record.adjoint!=(family==1))throw std::runtime_error("radial L0 source/adjoint flag disagrees with its synthesis set");
      if(record.rawStrain.empty()||(l0Seen_[key]&(1u<<family)))throw std::runtime_error("duplicate or empty radial L0 response record");
      const std::uint64_t offset=static_cast<std::uint64_t>(stream_.tellp());
      const auto &v=record.rawStrain.front();
      if(!finite(v[0])||!finite(v[1]))throw std::runtime_error("radial L0 response contains nonfinite values");
      writeValue(v[0]);writeValue(v[1]);
      writeIndex(l0OffsetsOffset_+(key*2+family)*sizeof(std::uint64_t),offset);
      l0Seen_[key]|=static_cast<std::uint8_t>(1u<<family);
    } else {
      const std::size_t family=record.adjoint?1:0;
      if((!record.adjoint&&record.synthesisSet!=0)||(record.adjoint&&record.synthesisSet!=1))
        throw std::runtime_error("unexpected synthesis set in response spool");
      const std::size_t index=key*2+family;
      if(recordSeen_[index])throw std::runtime_error("duplicate pole response record");
      if(record.labels.size()!=record.rawStrain.size()||record.labels.size()>static_cast<std::size_t>(5*(lmax_+1)))
        throw std::runtime_error("pole response labels and coefficients differ in length or exceed the degree bound");
      const std::uint64_t offset=static_cast<std::uint64_t>(stream_.tellp());
      writeValue(static_cast<std::uint64_t>(record.labels.size()));
      for(const auto &axis:record.frame.axes)for(double value:axis){if(!std::isfinite(value))throw std::runtime_error("pole record frame contains nonfinite values");writeValue(value);}
      for(std::size_t i=0;i<record.labels.size();++i){
        const auto [mp,degree]=record.labels[i];
        if(degree<0||degree>lmax_||mp<-2||mp>2||std::abs(mp)>degree||
           (i&&!(record.labels[i-1]<record.labels[i])))throw std::runtime_error("pole response labels are invalid or unordered");
        writeValue(static_cast<std::int32_t>(record.labels[i].first));writeValue(static_cast<std::int32_t>(record.labels[i].second));
        for(const Complex value:record.rawStrain[i]){if(!finite(value))throw std::runtime_error("pole response coefficient is nonfinite");writeValue(value);}
      }
      writeIndex(recordOffsetsOffset_+index*sizeof(std::uint64_t),offset);
      recordSeen_[index]=1;
    }
    const auto end=static_cast<std::uint64_t>(stream_.tellp());
    if(end>byteCap_)throw std::runtime_error("pole response spool exceeded configured byte cap");
  }

  void finalize() {
    requireOpen();
    for(auto seen:recordSeen_)if(!seen)throw std::runtime_error("response spool is missing a source or unit-adjoint record");
    for(auto seen:l0Seen_)if(seen!=3)throw std::runtime_error("response spool is missing a radial L0 record");
    for(std::size_t bin=0;bin<auditSet_.size();++bin) {
      if(!auditSet_[bin])throw std::runtime_error("response spool is missing a frequency adjoint audit");
      if(audits_[bin].frequencyBin!=static_cast<std::int64_t>(bin))
        throw std::runtime_error("response spool audits are not contiguous in the declared DC-inclusive range");
    }
    const auto resume=static_cast<std::uint64_t>(stream_.tellp());
    stream_.seekp(static_cast<std::streamoff>(auditOffset_));
    for(const auto &a:audits_){writeValue(a.frequencyBin);writeValue(a.z);writeValue(a.adjointMixing);writeValue(a.basisReconstructionRelative);}
    stream_.seekp(static_cast<std::streamoff>(resume));stream_.flush();
    if(!stream_)throw std::runtime_error("failed finalizing pole response spool");
    closed_=true;stream_.close();
  }

  std::uint64_t bytesWritten()const{return closed_?std::filesystem::file_size(path_):static_cast<std::uint64_t>(stream_.tellp());}
private:
  static bool finite(Complex x){return std::isfinite(x.real())&&std::isfinite(x.imag());}
  template<class T>void writeValue(const T &v){write(&v,sizeof(T));}
  void write(const void *p,std::size_t n){stream_.write(static_cast<const char*>(p),static_cast<std::streamsize>(n));if(!stream_)throw std::runtime_error("pole response spool write failed");}
  void writeIndex(std::uint64_t address,std::uint64_t value){const auto resume=static_cast<std::uint64_t>(stream_.tellp());stream_.seekp(static_cast<std::streamoff>(address));writeValue(value);stream_.seekp(static_cast<std::streamoff>(resume));if(!stream_)throw std::runtime_error("pole response index update failed");}
  void requireOpen()const{if(closed_||!stream_)throw std::runtime_error("pole response spool is closed");}
  std::string path_;int lmax_;std::size_t bins_;std::vector<NativeSide>sides_;std::uint64_t byteCap_;
  std::map<std::pair<int,std::uint64_t>,std::size_t>sideIndex_;std::vector<std::uint8_t>recordSeen_,l0Seen_;
  std::vector<SpoolAudit>audits_;std::vector<bool>auditSet_;std::fstream stream_;
  std::uint64_t recordOffsetsOffset_=0,l0OffsetsOffset_=0,auditOffset_=0,dataOffset_=0;bool closed_=false;
};

class PoleResponseSpoolReader {
public:
  explicit PoleResponseSpoolReader(const std::string &path):path_(path),stream_(path,std::ios::binary) {
    if(!stream_)throw std::runtime_error("cannot open pole response spool");
    fileBytes_=std::filesystem::file_size(path_);
    if(fileBytes_>static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()))
      throw std::runtime_error("pole response spool exceeds addressable stream offsets");
    std::array<char,8> magic{};read(magic.data(),magic.size());
    const std::array<char,8> expected{'D','S','G','S','H','S','P','1'};
    if(magic!=expected||readValue<std::uint32_t>()!=1)throw std::runtime_error("unsupported pole response spool version");
    lmax_=readValue<std::int32_t>();const std::uint64_t bins64=readValue<std::uint64_t>();const std::uint64_t sides64=readValue<std::uint64_t>();
    if(lmax_<1||lmax_>1536||bins64==0||bins64>static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max())+1||
       sides64==0||bins64>std::numeric_limits<std::size_t>::max()||sides64>std::numeric_limits<std::size_t>::max()||
       bins64>std::numeric_limits<std::size_t>::max()/sides64/2)
      throw std::runtime_error("invalid pole response spool header");
    const std::uint64_t keys=checkedMultiply(bins64,sides64,"spool key count overflows");
    const std::uint64_t indexBytes=checkedMultiply(checkedMultiply(keys,2,"spool index count overflows"),sizeof(std::uint64_t),"spool index bytes overflow");
    const std::uint64_t sideBytes=checkedMultiply(sides64,2*sizeof(std::int32_t)+3*sizeof(double)+sizeof(std::uint8_t),"native-side header bytes overflow");
    const std::uint64_t auditBytes=checkedMultiply(bins64,sizeof(std::int32_t)+2*sizeof(Complex)+sizeof(double),"audit table bytes overflow");
    if(sideBytes>512ull*1024ull*1024ull||auditBytes>512ull*1024ull*1024ull-sideBytes)
      throw std::runtime_error("response-spool side and audit metadata exceed their 512 MiB resident cap");
    std::uint64_t predictedDataOffset=checkedAdd(32,sideBytes,"spool header size overflows");
    predictedDataOffset=checkedAdd(predictedDataOffset,checkedMultiply(indexBytes,2,"spool index sections overflow"),"spool index end overflows");
    predictedDataOffset=checkedAdd(predictedDataOffset,auditBytes,"spool data offset overflows");
    if(predictedDataOffset>fileBytes_||predictedDataOffset>static_cast<std::uint64_t>(std::numeric_limits<std::streamoff>::max()))
      throw std::runtime_error("pole response spool file is truncated or header counts exceed its length");
    bins_=static_cast<std::size_t>(bins64);const std::size_t nSides=static_cast<std::size_t>(sides64);
    sides_.resize(nSides);for(auto &s:sides_){s.element=readValue<std::int32_t>();s.node=readValue<std::int32_t>();s.radiusMetres=readValue<double>();s.radiusModel=readValue<double>();s.widthMetres=readValue<double>();const auto fluid=readValue<std::uint8_t>();s.fluid=fluid!=0;
      if(s.element<0||s.node<0||fluid>1||!(std::isfinite(s.radiusMetres)&&s.radiusMetres>0.0&&
         std::isfinite(s.radiusModel)&&s.radiusModel>0.0&&std::isfinite(s.widthMetres)&&s.widthMetres>0.0))
        throw std::runtime_error("pole response spool native-side metadata is invalid");}
    recordOffsetsOffset_=static_cast<std::uint64_t>(stream_.tellg());
    recordOffsetsEnd_=checkedAdd(recordOffsetsOffset_,indexBytes,"record index end overflows");
    l0OffsetsOffset_=recordOffsetsEnd_;
    auditOffset_=checkedAdd(l0OffsetsOffset_,indexBytes,"L0 index end overflows");
    dataOffset_=checkedAdd(auditOffset_,auditBytes,"spool data offset overflows");
    if(dataOffset_!=predictedDataOffset)throw std::runtime_error("pole response spool header arithmetic is inconsistent");
    stream_.seekg(static_cast<std::streamoff>(auditOffset_));
    audits_.resize(static_cast<std::size_t>(bins_));for(auto &a:audits_){a.frequencyBin=readValue<std::int32_t>();a.z=readValue<Complex>();a.adjointMixing=readValue<Complex>();a.basisReconstructionRelative=readValue<double>();}
    for(std::size_t i=0;i<audits_.size();++i) {
      const auto &a=audits_[i];
      if(a.frequencyBin!=static_cast<std::int64_t>(i)||!finite(a.z)||!finite(a.adjointMixing)||
         !std::isfinite(a.basisReconstructionRelative)||a.basisReconstructionRelative<0.0)
        throw std::runtime_error("pole response spool audit identities or values are invalid");
    }
    validateAllOffsets(keys);
  }
  std::size_t frequencyCount()const noexcept{return static_cast<std::size_t>(bins_);}
  int lmax()const noexcept{return lmax_;}
  const std::vector<NativeSide>&sides()const noexcept{return sides_;}
  const SpoolAudit&frequencyAudit(std::size_t bin)const{return audits_.at(bin);}
  FrequencySideRecords read(std::size_t bin,std::size_t side) {
    if(bin>=bins_||side>=sides_.size())throw std::out_of_range("pole response spool key");
    const std::size_t key=bin*sides_.size()+side;
    FrequencySideRecords out;out.audit=audits_[bin];
    const std::uint64_t f0=readIndex(l0OffsetsOffset_+(key*2)*sizeof(std::uint64_t));
    const std::uint64_t f1=readIndex(l0OffsetsOffset_+(key*2+1)*sizeof(std::uint64_t));
    if(f0==0||f1==0)throw std::runtime_error("incomplete radial L0 spool records");
    seek(f0);out.l0.forwardRr=readValue<Complex>();out.l0.forwardHorizontalTrace=readValue<Complex>();
    seek(f1);out.l0.adjointRr=readValue<Complex>();out.l0.adjointHorizontalTrace=readValue<Complex>();
    if(!finite(out.l0.forwardRr)||!finite(out.l0.forwardHorizontalTrace)||
       !finite(out.l0.adjointRr)||!finite(out.l0.adjointHorizontalTrace))
      throw std::runtime_error("radial L0 spool values are nonfinite");
    const auto o0=readIndex(recordOffsetsOffset_+(key*2)*sizeof(std::uint64_t));
    const auto o1=readIndex(recordOffsetsOffset_+(key*2+1)*sizeof(std::uint64_t));
    if(o0==0||o1==0)throw std::runtime_error("incomplete pole response spool");
    readRecord(o0,side,false,out.forward);readRecord(o1,side,true,out.unitAdjoint);
    return out;
  }
private:
  struct Interval {std::uint64_t begin,end;};
  static std::uint64_t checkedMultiply(std::uint64_t a,std::uint64_t b,const char *what) {
    if(b&&a>std::numeric_limits<std::uint64_t>::max()/b)throw std::runtime_error(what);return a*b;
  }
  static std::uint64_t checkedAdd(std::uint64_t a,std::uint64_t b,const char *what) {
    if(a>std::numeric_limits<std::uint64_t>::max()-b)throw std::runtime_error(what);return a+b;
  }
  static bool finite(Complex z){return std::isfinite(z.real())&&std::isfinite(z.imag());}
  template<class T>T readValue(){T v{};read(&v,sizeof(T));return v;}
  void read(void *p,std::size_t n){stream_.read(static_cast<char*>(p),static_cast<std::streamsize>(n));if(!stream_)throw std::runtime_error("pole response spool read failed");}
  void seek(std::uint64_t offset){stream_.clear();stream_.seekg(static_cast<std::streamoff>(offset));if(!stream_)throw std::runtime_error("pole response spool seek failed");}
  std::uint64_t readIndex(std::uint64_t address) {
    if(address>fileBytes_-sizeof(std::uint64_t)||
       !((address>=recordOffsetsOffset_&&address+sizeof(std::uint64_t)<=recordOffsetsEnd_)||
         (address>=l0OffsetsOffset_&&address+sizeof(std::uint64_t)<=auditOffset_)))
      throw std::runtime_error("pole response index-table address is outside its declared region");
    seek(address);return readValue<std::uint64_t>();
  }
  void validateAllOffsets(std::uint64_t keys) {
    constexpr std::uint64_t maxValidationMemory=256ull*1024ull*1024ull;
    const std::uint64_t count=checkedMultiply(keys,4,"spool interval count overflows");
    const std::uint64_t bytes=checkedMultiply(count,sizeof(Interval),"spool interval memory overflows");
    if(bytes>maxValidationMemory||count>std::numeric_limits<std::size_t>::max())
      throw std::runtime_error("spool offset validation exceeds its 256 MiB temporary memory cap");
    std::vector<Interval> intervals;intervals.reserve(static_cast<std::size_t>(count));
    for(std::uint64_t key=0;key<keys;++key) {
      for(std::size_t family=0;family<2;++family) {
        const auto address=recordOffsetsOffset_+(key*2+family)*sizeof(std::uint64_t);
        const auto offset=readIndex(address);if(offset==0)throw std::runtime_error("response spool has an empty pole offset");
        if(offset<dataOffset_||offset>fileBytes_-sizeof(std::uint64_t))throw std::runtime_error("pole response offset is outside the data region");
        seek(offset);const auto n=readValue<std::uint64_t>();
        if(n>static_cast<std::uint64_t>(5*(lmax_+1)))throw std::runtime_error("pole response label count exceeds its degree bound");
        const auto labelBytes=checkedMultiply(n,2*sizeof(std::int32_t)+6*sizeof(Complex),"pole record length overflows");
        const auto length=checkedAdd(sizeof(std::uint64_t)+9*sizeof(double),labelBytes,"pole record length overflows");
        const auto end=checkedAdd(offset,length,"pole response end overflows");
        if(end>fileBytes_)throw std::runtime_error("pole response record extends past the spool file");
        intervals.push_back({offset,end});
      }
      for(std::size_t family=0;family<2;++family) {
        const auto address=l0OffsetsOffset_+(key*2+family)*sizeof(std::uint64_t);
        const auto offset=readIndex(address);
        if(offset<dataOffset_||offset>fileBytes_||2*sizeof(Complex)>fileBytes_-offset)
          throw std::runtime_error("radial L0 offset or payload is outside the spool data region");
        intervals.push_back({offset,offset+2*sizeof(Complex)});
      }
    }
    std::sort(intervals.begin(),intervals.end(),[](const Interval &a,const Interval &b){return a.begin<b.begin;});
    for(std::size_t i=1;i<intervals.size();++i)
      if(intervals[i].begin<intervals[i-1].end)
        throw std::runtime_error("pole response spool payload offsets overlap or alias another record");
  }
  void readRecord(std::uint64_t offset,std::size_t side,bool adjoint,PoleRecord &record) {
    seek(offset);const std::uint64_t n=readValue<std::uint64_t>();
    if(n>static_cast<std::uint64_t>(5*(lmax_+1)))throw std::runtime_error("pole response record label count is invalid");
    record.adjoint=adjoint;record.labels.resize(static_cast<std::size_t>(n));record.rawStrain.resize(static_cast<std::size_t>(n));
    for(auto &axis:record.frame.axes)for(double &value:axis){value=readValue<double>();if(!std::isfinite(value))throw std::runtime_error("pole record frame contains nonfinite data");}
    for(std::size_t i=0;i<n;++i){
      record.labels[i].first=readValue<std::int32_t>();record.labels[i].second=readValue<std::int32_t>();
      const auto [mp,degree]=record.labels[i];
      if(degree<0||degree>lmax_||mp<-2||mp>2||std::abs(mp)>degree||
         (i&&!(record.labels[i-1]<record.labels[i])))throw std::runtime_error("pole response labels are invalid or unordered");
      for(auto &v:record.rawStrain[i]){v=readValue<Complex>();if(!finite(v))throw std::runtime_error("pole response coefficient is nonfinite");}
    }
    record.element=sides_[side].element;record.radiusModel=sides_[side].radiusModel;
    record.synthesisSet=adjoint?1:0;
  }
  std::string path_;std::ifstream stream_;int lmax_=0;std::size_t bins_=0;std::uint64_t fileBytes_=0;
  std::uint64_t recordOffsetsOffset_=0,recordOffsetsEnd_=0,l0OffsetsOffset_=0,auditOffset_=0,dataOffset_=0;
  std::vector<NativeSide>sides_;std::vector<SpoolAudit>audits_;
};
} // namespace dspec_global_gsh
