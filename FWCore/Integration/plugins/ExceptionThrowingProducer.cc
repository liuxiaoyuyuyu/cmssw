// -*- C++ -*-
//
// Package:     FWCore/Integration
// Class  :     ExceptionThrowingProducer
//
// Implementation:
//     Intentionally throws exceptions in various Framework transitions.
//     You can configure which transition. Includes some tests of
//     the Framework behavior after an exception occurs.
//
// Original Author:  W. David Dagenhart
//         Created:  26 September 2022

#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Integration/plugins/TestServiceOne.h"
#include "FWCore/Integration/plugins/TestServiceTwo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include <atomic>
#include <limits>
#include <memory>

constexpr unsigned int kTestStreams = 4;
constexpr unsigned int kUnset = std::numeric_limits<unsigned int>::max();
constexpr unsigned int kNumberOfTestModules = 2;
// 3 modules configured plus TriggerResultsInserter plus PathStatusInserter plus EndPathStatusInserter
constexpr unsigned int kNumberOfModulesWithJobTransitions = 6;
constexpr unsigned int kNumberOfModulesWithStreamTransitions = kNumberOfModulesWithJobTransitions;
constexpr unsigned int kNumberOfModulesWithProcessBlockTransitions = 2;
constexpr unsigned int kStreamToThrowOn = 2;

namespace edmtest {

  namespace {
    struct Cache {};
  }  // namespace

  class ExceptionThrowingProducer : public edm::global::EDProducer<edm::StreamCache<Cache>,
                                                                   edm::RunCache<Cache>,
                                                                   edm::LuminosityBlockCache<Cache>,
                                                                   edm::WatchProcessBlock> {
  public:
    explicit ExceptionThrowingProducer(edm::ParameterSet const&);

    ~ExceptionThrowingProducer() noexcept(false) override;
    ExceptionThrowingProducer(const ExceptionThrowingProducer&) = delete;
    ExceptionThrowingProducer& operator=(const ExceptionThrowingProducer&) = delete;
    ExceptionThrowingProducer(ExceptionThrowingProducer&&) = delete;
    ExceptionThrowingProducer& operator=(ExceptionThrowingProducer&&) = delete;

    void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

    std::shared_ptr<Cache> globalBeginRun(edm::Run const&, edm::EventSetup const&) const override;
    void globalEndRun(edm::Run const&, edm::EventSetup const&) const override;
    std::shared_ptr<Cache> globalBeginLuminosityBlock(edm::LuminosityBlock const&,
                                                      edm::EventSetup const&) const override;
    void globalEndLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) const override;

    std::unique_ptr<Cache> beginStream(edm::StreamID) const override;
    void streamBeginRun(edm::StreamID, edm::Run const&, edm::EventSetup const&) const override;
    void streamBeginLuminosityBlock(edm::StreamID, edm::LuminosityBlock const&, edm::EventSetup const&) const override;
    void streamEndLuminosityBlock(edm::StreamID, edm::LuminosityBlock const&, edm::EventSetup const&) const override;
    void streamEndRun(edm::StreamID, edm::Run const&, edm::EventSetup const&) const override;
    void endStream(edm::StreamID) const override;

    void beginProcessBlock(edm::ProcessBlock const& processBlock) override;
    void endProcessBlock(edm::ProcessBlock const& processBlock) override;

    void beginJob() override;
    void endJob() override;

    static void fillDescriptions(edm::ConfigurationDescriptions&);

  private:
    const bool verbose_;

    const edm::EventID eventIDThrowOnEvent_;
    const edm::EventID eventIDThrowOnGlobalBeginRun_;
    const edm::EventID eventIDThrowOnGlobalBeginLumi_;
    const edm::EventID eventIDThrowOnGlobalEndRun_;
    const edm::EventID eventIDThrowOnGlobalEndLumi_;
    const edm::EventID eventIDThrowOnStreamBeginRun_;
    const edm::EventID eventIDThrowOnStreamBeginLumi_;
    const edm::EventID eventIDThrowOnStreamEndRun_;
    const edm::EventID eventIDThrowOnStreamEndLumi_;
    const bool throwInBeginJob_;
    const bool throwInBeginStream_;
    const bool throwInBeginProcessBlock_;
    const bool throwInEndProcessBlock_;
    const bool throwInEndStream_;
    const bool throwInEndJob_;

    unsigned int nBeginJob_ = 0;
    mutable std::atomic<unsigned int> nBeginStream_{0};
    unsigned int nBeginProcessBlock_ = 0;
    unsigned int nEndProcessBlock_ = 0;
    mutable std::atomic<unsigned int> nEndStream_{0};
    unsigned int nEndJob_ = 0;

    mutable std::vector<unsigned int> nStreamBeginLumi_;
    mutable std::vector<unsigned int> nStreamEndLumi_;
    mutable std::atomic<unsigned int> nGlobalBeginLumi_{0};
    mutable std::atomic<unsigned int> nGlobalEndLumi_{0};

    mutable std::vector<unsigned int> nStreamBeginRun_;
    mutable std::vector<unsigned int> nStreamEndRun_;
    mutable std::atomic<unsigned int> nGlobalBeginRun_{0};
    mutable std::atomic<unsigned int> nGlobalEndRun_{0};

    const unsigned int expectedNBeginJob_;
    const unsigned int expectedNBeginStream_;
    const unsigned int expectedNBeginProcessBlock_;
    const unsigned int expectedNEndProcessBlock_;
    const unsigned int expectedNEndStream_;
    const unsigned int expectedNEndJob_;
    const bool expectNoRunsProcessed_;

    const unsigned int expectedOffsetNoEndJob_;
    const unsigned int expectedOffsetNoEndStream_;
    const unsigned int expectedOffsetNoEndProcessBlock_;

    const unsigned int expectedStreamBeginLumi_;
    const unsigned int expectedOffsetNoStreamEndLumi_;
    mutable unsigned int streamWithBeginLumiException_ = kUnset;
    const unsigned int expectedGlobalBeginLumi_;
    const unsigned int expectedOffsetNoGlobalEndLumi_;
    const unsigned int expectedOffsetNoWriteLumi_;

    const unsigned int expectedStreamBeginRun_;
    const unsigned int expectedOffsetNoStreamEndRun_;
    mutable unsigned int streamWithBeginRunException_ = kUnset;
    const unsigned int expectedGlobalBeginRun_;
    const unsigned int expectedOffsetNoGlobalEndRun_;
    const unsigned int expectedOffsetNoWriteRun_;

    mutable std::atomic<bool> streamBeginLumiExceptionOccurred_ = false;
    mutable std::atomic<bool> streamEndLumiExceptionOccurred_ = false;
    mutable std::atomic<bool> globalBeginLumiExceptionOccurred_ = false;

    mutable std::atomic<bool> streamBeginRunExceptionOccurred_ = false;
    mutable std::atomic<bool> streamEndRunExceptionOccurred_ = false;
    mutable std::atomic<bool> globalBeginRunExceptionOccurred_ = false;
  };

  ExceptionThrowingProducer::ExceptionThrowingProducer(edm::ParameterSet const& pset)
      : verbose_(pset.getUntrackedParameter<bool>("verbose")),
        eventIDThrowOnEvent_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnEvent")),
        eventIDThrowOnGlobalBeginRun_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnGlobalBeginRun")),
        eventIDThrowOnGlobalBeginLumi_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnGlobalBeginLumi")),
        eventIDThrowOnGlobalEndRun_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnGlobalEndRun")),
        eventIDThrowOnGlobalEndLumi_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnGlobalEndLumi")),
        eventIDThrowOnStreamBeginRun_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnStreamBeginRun")),
        eventIDThrowOnStreamBeginLumi_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnStreamBeginLumi")),
        eventIDThrowOnStreamEndRun_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnStreamEndRun")),
        eventIDThrowOnStreamEndLumi_(pset.getUntrackedParameter<edm::EventID>("eventIDThrowOnStreamEndLumi")),
        throwInBeginJob_(pset.getUntrackedParameter<bool>("throwInBeginJob")),
        throwInBeginStream_(pset.getUntrackedParameter<bool>("throwInBeginStream")),
        throwInBeginProcessBlock_(pset.getUntrackedParameter<bool>("throwInBeginProcessBlock")),
        throwInEndProcessBlock_(pset.getUntrackedParameter<bool>("throwInEndProcessBlock")),
        throwInEndStream_(pset.getUntrackedParameter<bool>("throwInEndStream")),
        throwInEndJob_(pset.getUntrackedParameter<bool>("throwInEndJob")),
        nStreamBeginLumi_(kTestStreams, 0),
        nStreamEndLumi_(kTestStreams, 0),
        nStreamBeginRun_(kTestStreams, 0),
        nStreamEndRun_(kTestStreams, 0),
        expectedNBeginJob_(pset.getUntrackedParameter<unsigned int>("expectedNBeginJob")),
        expectedNBeginStream_(pset.getUntrackedParameter<unsigned int>("expectedNBeginStream")),
        expectedNBeginProcessBlock_(pset.getUntrackedParameter<unsigned int>("expectedNBeginProcessBlock")),
        expectedNEndProcessBlock_(pset.getUntrackedParameter<unsigned int>("expectedNEndProcessBlock")),
        expectedNEndStream_(pset.getUntrackedParameter<unsigned int>("expectedNEndStream")),
        expectedNEndJob_(pset.getUntrackedParameter<unsigned int>("expectedNEndJob")),
        expectNoRunsProcessed_(pset.getUntrackedParameter<bool>("expectNoRunsProcessed")),
        expectedOffsetNoEndJob_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoEndJob")),
        expectedOffsetNoEndStream_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoEndStream")),
        expectedOffsetNoEndProcessBlock_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoEndProcessBlock")),
        expectedStreamBeginLumi_(pset.getUntrackedParameter<unsigned int>("expectedStreamBeginLumi")),
        expectedOffsetNoStreamEndLumi_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoStreamEndLumi")),
        expectedGlobalBeginLumi_(pset.getUntrackedParameter<unsigned int>("expectedGlobalBeginLumi")),
        expectedOffsetNoGlobalEndLumi_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoGlobalEndLumi")),
        expectedOffsetNoWriteLumi_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoWriteLumi")),
        expectedStreamBeginRun_(pset.getUntrackedParameter<unsigned int>("expectedStreamBeginRun")),
        expectedOffsetNoStreamEndRun_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoStreamEndRun")),
        expectedGlobalBeginRun_(pset.getUntrackedParameter<unsigned int>("expectedGlobalBeginRun")),
        expectedOffsetNoGlobalEndRun_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoGlobalEndRun")),
        expectedOffsetNoWriteRun_(pset.getUntrackedParameter<unsigned int>("expectedOffsetNoWriteRun")) {}

  ExceptionThrowingProducer::~ExceptionThrowingProducer() noexcept(false) {
    if (nBeginJob_ != expectedNBeginJob_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "ExceptionThrowingProducer::beginJob did not run expected number of times.";
      edm::LogAbsolute("ExceptionThrowingProducer") << "At least one test in ExceptionThrowingProducer FAILED";
    } else if (nEndJob_ != expectedNEndJob_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "ExceptionThrowingProducer::endJob did not run expected number of times.";
      edm::LogAbsolute("ExceptionThrowingProducer") << "At least one test in ExceptionThrowingProducer FAILED";
    } else if (nEndJob_ == 0u && nBeginProcessBlock_ == 0 && nEndProcessBlock_ == 0 && nBeginStream_.load() == 0 &&
               nEndStream_.load() == 0 && nGlobalBeginRun_.load() == 0 && nGlobalEndRun_.load() == 0) {
      edm::LogAbsolute("ExceptionThrowingProducer") << "All tests in ExceptionThrowingProducer PASSED";
    }

    edm::Service<edmtest::TestServiceOne> serviceOne;

    if (verbose_) {
      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPostEndJob = " << serviceOne->nPostEndJob();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleEndJob = " << serviceOne->nPreModuleEndJob();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleEndJob = " << serviceOne->nPostModuleEndJob();
    }

    if (serviceOne->nPostEndJob() != nBeginJob_ ||
        serviceOne->nPreModuleEndJob() != nBeginJob_ * kNumberOfModulesWithJobTransitions - expectedOffsetNoEndJob_ ||
        serviceOne->nPostModuleEndJob() != nBeginJob_ * kNumberOfModulesWithJobTransitions - expectedOffsetNoEndJob_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service end job related signals in TestServiceOne, checking while in "
             "ExceptionThrowingProducer destructor";
      edm::LogAbsolute("ExceptionThrowingProducer") << "At least one test in ExceptionThrowingProducer FAILED";
    }
  }

  void ExceptionThrowingProducer::produce(edm::StreamID, edm::Event& event, edm::EventSetup const&) const {
    if (event.id() == eventIDThrowOnEvent_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::produce, module configured to throw on: " << eventIDThrowOnEvent_;
    }
  }

  std::shared_ptr<Cache> ExceptionThrowingProducer::globalBeginRun(edm::Run const& run, edm::EventSetup const&) const {
    ++nGlobalBeginRun_;
    if (edm::EventID(run.id().run(), edm::invalidLuminosityBlockNumber, edm::invalidEventNumber) ==
        eventIDThrowOnGlobalBeginRun_) {
      globalBeginRunExceptionOccurred_.store(true);
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::globalBeginRun, module configured to throw on: "
          << eventIDThrowOnGlobalBeginRun_;
    }
    return std::make_shared<Cache>();
  }

  void ExceptionThrowingProducer::globalEndRun(edm::Run const& run, edm::EventSetup const&) const {
    ++nGlobalEndRun_;
    if (edm::EventID(run.id().run(), edm::invalidLuminosityBlockNumber, edm::invalidEventNumber) ==
        eventIDThrowOnGlobalEndRun_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::globalEndRun, module configured to throw on: " << eventIDThrowOnGlobalEndRun_;
    }
  }

  std::shared_ptr<Cache> ExceptionThrowingProducer::globalBeginLuminosityBlock(edm::LuminosityBlock const& lumi,
                                                                               edm::EventSetup const&) const {
    ++nGlobalBeginLumi_;
    if (edm::EventID(lumi.id().run(), lumi.id().luminosityBlock(), edm::invalidEventNumber) ==
        eventIDThrowOnGlobalBeginLumi_) {
      globalBeginLumiExceptionOccurred_.store(true);
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::globalBeginLuminosityBlock, module configured to throw on: "
          << eventIDThrowOnGlobalBeginLumi_;
    }
    return std::make_shared<Cache>();
  }

  void ExceptionThrowingProducer::globalEndLuminosityBlock(edm::LuminosityBlock const& lumi,
                                                           edm::EventSetup const&) const {
    ++nGlobalEndLumi_;
    if (edm::EventID(lumi.id().run(), lumi.id().luminosityBlock(), edm::invalidEventNumber) ==
        eventIDThrowOnGlobalEndLumi_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::globalEndLuminosityBlock, module configured to throw on: "
          << eventIDThrowOnGlobalEndLumi_;
    }
  }

  std::unique_ptr<Cache> ExceptionThrowingProducer::beginStream(edm::StreamID streamID) const {
    ++nBeginStream_;
    if (throwInBeginStream_ && streamID.value() == kStreamToThrowOn) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::beginStream, module configured to throw during beginStream for stream: "
          << kStreamToThrowOn;
    }
    return std::make_unique<Cache>();
  }

  void ExceptionThrowingProducer::streamBeginRun(edm::StreamID iStream,
                                                 edm::Run const& run,
                                                 edm::EventSetup const&) const {
    if (iStream < kTestStreams) {
      ++nStreamBeginRun_[iStream];
    }

    bool expected = false;
    if (edm::EventID(run.id().run(), edm::invalidLuminosityBlockNumber, edm::invalidEventNumber) ==
            eventIDThrowOnStreamBeginRun_ &&
        streamBeginRunExceptionOccurred_.compare_exchange_strong(expected, true)) {
      // Remember which stream threw
      streamWithBeginRunException_ = iStream.value();
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::streamBeginRun, module configured to throw on: "
          << eventIDThrowOnStreamBeginRun_;
    }
  }

  void ExceptionThrowingProducer::streamBeginLuminosityBlock(edm::StreamID iStream,
                                                             edm::LuminosityBlock const& lumi,
                                                             edm::EventSetup const&) const {
    if (iStream < kTestStreams) {
      ++nStreamBeginLumi_[iStream];
    }

    // Throw if this lumi's ID matches the configured ID (this code is written so
    // that only the first stream to match it will throw).
    bool expected = false;
    if (edm::EventID(lumi.run(), lumi.id().luminosityBlock(), edm::invalidEventNumber) ==
            eventIDThrowOnStreamBeginLumi_ &&
        streamBeginLumiExceptionOccurred_.compare_exchange_strong(expected, true)) {
      // Remember which stream threw
      streamWithBeginLumiException_ = iStream.value();

      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::streamBeginLuminosityBlock, module configured to throw on: "
          << eventIDThrowOnStreamBeginLumi_;
    }
  }

  void ExceptionThrowingProducer::streamEndLuminosityBlock(edm::StreamID iStream,
                                                           edm::LuminosityBlock const& lumi,
                                                           edm::EventSetup const&) const {
    if (iStream < kTestStreams) {
      ++nStreamEndLumi_[iStream];
    }

    bool expected = false;
    if (edm::EventID(lumi.run(), lumi.id().luminosityBlock(), edm::invalidEventNumber) ==
            eventIDThrowOnStreamEndLumi_ &&
        streamEndLumiExceptionOccurred_.compare_exchange_strong(expected, true)) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::streamEndLuminosityBlock, module configured to throw on: "
          << eventIDThrowOnStreamEndLumi_;
    }
  }

  void ExceptionThrowingProducer::streamEndRun(edm::StreamID iStream,
                                               edm::Run const& run,
                                               edm::EventSetup const&) const {
    if (iStream < kTestStreams) {
      ++nStreamEndRun_[iStream];
    }

    bool expected = false;
    if (edm::EventID(run.id().run(), edm::invalidLuminosityBlockNumber, edm::invalidEventNumber) ==
            eventIDThrowOnStreamEndRun_ &&
        streamEndRunExceptionOccurred_.compare_exchange_strong(expected, true)) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::streamEndRun, module configured to throw on: " << eventIDThrowOnStreamEndRun_;
    }
  }

  void ExceptionThrowingProducer::endStream(edm::StreamID streamID) const {
    ++nEndStream_;
    if (throwInEndStream_ && streamID.value() == kStreamToThrowOn) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::endStream, module configured to throw during endStream for stream: "
          << kStreamToThrowOn;
    }
  }

  void ExceptionThrowingProducer::beginProcessBlock(edm::ProcessBlock const& processBlock) {
    ++nBeginProcessBlock_;
    if (throwInBeginProcessBlock_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::beginProcessBlock, module configured to throw during beginProcessBlock";
    }
  }

  void ExceptionThrowingProducer::endProcessBlock(edm::ProcessBlock const& processBlock) {
    ++nEndProcessBlock_;
    if (throwInEndProcessBlock_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::endProcessBlock, module configured to throw during endProcessBlock";
    }
  }

  void ExceptionThrowingProducer::beginJob() {
    ++nBeginJob_;
    if (throwInBeginJob_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::beginJob, module configured to throw during beginJob";
    }
  }

  void ExceptionThrowingProducer::endJob() {
    ++nEndJob_;

    bool testsPass = true;

    if (nBeginStream_.load() != expectedNBeginStream_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of beginStream transitions, saw: " << nBeginStream_.load()
          << " expected: " << expectedNBeginStream_;
      testsPass = false;
    }

    if (nBeginProcessBlock_ != expectedNBeginProcessBlock_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of beginProcessBlock transitions, saw: " << nBeginProcessBlock_
          << " expected: " << expectedNBeginProcessBlock_;
      testsPass = false;
    }

    if (nEndProcessBlock_ != expectedNEndProcessBlock_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of endProcessBlock transitions, saw: " << nEndProcessBlock_
          << " expected: " << expectedNEndProcessBlock_;
      testsPass = false;
    }

    if (nEndStream_.load() != expectedNEndStream_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of endStream transitions, saw: " << nEndStream_.load()
          << " expected: " << expectedNEndStream_;
      testsPass = false;
    }

    if (expectNoRunsProcessed_ && nGlobalBeginRun_.load() != 0) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of globalBeginRun transitions, saw: " << nGlobalBeginRun_.load()
          << " expected: " << 0;
      testsPass = false;
    }

    unsigned int totalStreamBeginLumi = 0;
    unsigned int i = 0;
    for (auto const& nStreamBeginLumi : nStreamBeginLumi_) {
      totalStreamBeginLumi += nStreamBeginLumi;

      // Don't know exact number to expect because streams might skip a lumi so
      // only throw if it is greater than the maximum possible and we only know
      // that for sure if the exception was thrown in stream begin lumi.
      if (nStreamBeginLumi > expectedStreamBeginLumi_ && streamWithBeginLumiException_ != kUnset) {
        edm::LogAbsolute("ExceptionThrowingProducer")
            << "FAILED: More than maximum possible number of streamBeginLumi transitions, stream " << i << " saw "
            << nStreamBeginLumi << " max possible " << expectedStreamBeginLumi_;
        testsPass = false;
      }
      unsigned int expectedStreamEndLumi =
          (streamWithBeginLumiException_ == i) ? nStreamBeginLumi - 1 : nStreamBeginLumi;
      if (nStreamEndLumi_[i] != expectedStreamEndLumi) {
        edm::LogAbsolute("ExceptionThrowingProducer")
            << "FAILED: Unexpected number of streamEndLumi transitions, stream " << i << " saw " << nStreamEndLumi_[i]
            << " expected " << expectedStreamEndLumi;
        testsPass = false;
      }

      ++i;
    }

    unsigned int totalStreamBeginRun = 0;
    i = 0;
    for (auto const& nStreamBeginRun : nStreamBeginRun_) {
      totalStreamBeginRun += nStreamBeginRun;

      // Don't know exact number to expect because streams might skip a run (not yet
      // but probably in the future) so only throw if it is greater than the maximum
      // possible and we only know that for sure if the exception was thrown in stream begin run.
      if (nStreamBeginRun > expectedStreamBeginRun_ && streamWithBeginRunException_ != kUnset) {
        edm::LogAbsolute("ExceptionThrowingProducer")
            << "FAILED: More than maximum possible number of streamBeginRun transitions, stream " << i << " saw "
            << nStreamBeginRun << " max possible " << expectedStreamBeginRun_;
        testsPass = false;
      }
      unsigned int expectedStreamEndRun = (streamWithBeginRunException_ == i) ? nStreamBeginRun - 1 : nStreamBeginRun;
      if (nStreamEndRun_[i] != expectedStreamEndRun) {
        edm::LogAbsolute("ExceptionThrowingProducer")
            << "FAILED: Unexpected number of streamEndRun transitions, stream " << i << " saw " << nStreamEndRun_[i]
            << " expected " << expectedStreamEndRun;
        testsPass = false;
      }

      ++i;
    }

    // There has to be at least as many global begin lumi transitions
    // as expected. Because of concurrency, the Framework might already have
    // started other lumis ahead of the one where an exception occurs.
    if (expectedGlobalBeginLumi_ > 0 && nGlobalBeginLumi_.load() < expectedGlobalBeginLumi_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Less than the expected number of globalBeginLumi transitions, expected at least "
          << expectedGlobalBeginLumi_ << " saw " << nGlobalBeginLumi_.load();
      testsPass = false;
    }

    // There has to be at least as many global begin run transitions
    // as expected. Because of concurrency, the Framework might already have
    // started other runs ahead of the one where an exception occurs.
    if (expectedGlobalBeginRun_ > 0 && nGlobalBeginRun_.load() < expectedGlobalBeginRun_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Less than the expected number of globalBeginRun transitions, expected at least "
          << expectedGlobalBeginRun_ << " saw " << nGlobalBeginRun_.load();
      testsPass = false;
    }

    unsigned int expectedGlobalEndLumi =
        globalBeginLumiExceptionOccurred_.load() ? nGlobalBeginLumi_.load() - 1 : nGlobalBeginLumi_.load();
    if (nGlobalEndLumi_.load() != expectedGlobalEndLumi) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: number of global end lumi transitions not equal to expected value, expected "
          << expectedGlobalEndLumi << " saw " << nGlobalEndLumi_.load();
      testsPass = false;
    }

    unsigned int expectedGlobalEndRun =
        globalBeginRunExceptionOccurred_.load() ? nGlobalBeginRun_.load() - 1 : nGlobalBeginRun_.load();
    if (nGlobalEndRun_.load() != expectedGlobalEndRun) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: number of global end run transitions not equal to expected value, expected "
          << expectedGlobalEndRun << " saw " << nGlobalEndRun_.load();
      testsPass = false;
    }

    edm::Service<edmtest::TestServiceOne> serviceOne;
    if (serviceOne->nPreBeginJob() != nBeginJob_ || serviceOne->nPostBeginJob() != nBeginJob_ ||
        serviceOne->nPreEndJob() != nBeginJob_ ||
        serviceOne->nPreModuleBeginJob() != nBeginJob_ * kNumberOfModulesWithJobTransitions ||
        serviceOne->nPostModuleBeginJob() != nBeginJob_ * kNumberOfModulesWithJobTransitions) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, beginJob and endJob";
      testsPass = false;
    }

    if (serviceOne->nPreBeginStream() != nBeginStream_ || serviceOne->nPostBeginStream() != nBeginStream_ ||
        serviceOne->nPreEndStream() != nBeginStream_ || serviceOne->nPostEndStream() != nBeginStream_ ||
        serviceOne->nPreModuleBeginStream() != nBeginStream_ * kNumberOfModulesWithStreamTransitions ||
        serviceOne->nPostModuleBeginStream() != nBeginStream_ * kNumberOfModulesWithStreamTransitions ||
        serviceOne->nPreModuleEndStream() !=
            nBeginStream_ * kNumberOfModulesWithStreamTransitions - expectedOffsetNoEndStream_ ||
        serviceOne->nPostModuleEndStream() !=
            nBeginStream_ * kNumberOfModulesWithStreamTransitions - expectedOffsetNoEndStream_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, beginStream and endStream";
      testsPass = false;
    }

    if (serviceOne->nPreBeginProcessBlock() != nBeginProcessBlock_ ||
        serviceOne->nPostBeginProcessBlock() != nBeginProcessBlock_ ||
        serviceOne->nPreEndProcessBlock() != nBeginProcessBlock_ ||
        serviceOne->nPostEndProcessBlock() != nBeginProcessBlock_ ||
        serviceOne->nPreModuleBeginProcessBlock() !=
            nBeginProcessBlock_ * kNumberOfModulesWithProcessBlockTransitions ||
        serviceOne->nPostModuleBeginProcessBlock() !=
            nBeginProcessBlock_ * kNumberOfModulesWithProcessBlockTransitions ||
        serviceOne->nPreModuleEndProcessBlock() !=
            nBeginProcessBlock_ * kNumberOfModulesWithProcessBlockTransitions - expectedOffsetNoEndProcessBlock_ ||
        serviceOne->nPostModuleEndProcessBlock() !=
            nBeginProcessBlock_ * kNumberOfModulesWithProcessBlockTransitions - expectedOffsetNoEndProcessBlock_) {
      edm::LogAbsolute("ExceptionThrowingProducer") << "FAILED: Unexpected number of service transitions in "
                                                       "TestServiceOne, beginProcessBlock and endProcessBlock";
      testsPass = false;
    }

    if (serviceOne->nPreStreamBeginLumi() != totalStreamBeginLumi ||
        serviceOne->nPostStreamBeginLumi() != totalStreamBeginLumi ||
        serviceOne->nPreStreamEndLumi() != totalStreamBeginLumi ||
        serviceOne->nPostStreamEndLumi() != totalStreamBeginLumi ||
        serviceOne->nPreModuleStreamBeginLumi() != totalStreamBeginLumi * kNumberOfTestModules ||
        serviceOne->nPostModuleStreamBeginLumi() != totalStreamBeginLumi * kNumberOfTestModules ||
        serviceOne->nPreModuleStreamEndLumi() !=
            totalStreamBeginLumi * kNumberOfTestModules - expectedOffsetNoStreamEndLumi_ ||
        serviceOne->nPostModuleStreamEndLumi() !=
            totalStreamBeginLumi * kNumberOfTestModules - expectedOffsetNoStreamEndLumi_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, stream lumi";
      testsPass = false;
    }

    edm::Service<edmtest::TestServiceTwo> serviceTwo;
    if (serviceTwo->nPreStreamBeginLumi() != totalStreamBeginLumi ||
        serviceTwo->nPostStreamBeginLumi() != totalStreamBeginLumi ||
        serviceTwo->nPreStreamEndLumi() != totalStreamBeginLumi ||
        serviceTwo->nPostStreamEndLumi() != totalStreamBeginLumi ||
        serviceTwo->nPreModuleStreamBeginLumi() != totalStreamBeginLumi * kNumberOfTestModules ||
        serviceTwo->nPostModuleStreamBeginLumi() != totalStreamBeginLumi * kNumberOfTestModules ||
        serviceTwo->nPreModuleStreamEndLumi() !=
            totalStreamBeginLumi * kNumberOfTestModules - expectedOffsetNoStreamEndLumi_ ||
        serviceTwo->nPostModuleStreamEndLumi() !=
            totalStreamBeginLumi * kNumberOfTestModules - expectedOffsetNoStreamEndLumi_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceTwo, stream lumi";
      testsPass = false;
    }

    unsigned int nGlobalBeginLumi = nGlobalBeginLumi_.load();

    if (serviceOne->nPreGlobalBeginLumi() != nGlobalBeginLumi ||
        serviceOne->nPostGlobalBeginLumi() != nGlobalBeginLumi || serviceOne->nPreGlobalEndLumi() != nGlobalBeginLumi ||
        serviceOne->nPostGlobalEndLumi() != nGlobalBeginLumi ||
        serviceOne->nPreModuleGlobalBeginLumi() != nGlobalBeginLumi * kNumberOfTestModules ||
        serviceOne->nPostModuleGlobalBeginLumi() != nGlobalBeginLumi * kNumberOfTestModules ||
        serviceOne->nPreModuleGlobalEndLumi() !=
            nGlobalBeginLumi * kNumberOfTestModules - expectedOffsetNoGlobalEndLumi_ ||
        serviceOne->nPostModuleGlobalEndLumi() !=
            nGlobalBeginLumi * kNumberOfTestModules - expectedOffsetNoGlobalEndLumi_ ||
        serviceOne->nPreGlobalWriteLumi() != nGlobalBeginLumi - expectedOffsetNoWriteLumi_ ||
        serviceOne->nPostGlobalWriteLumi() != nGlobalBeginLumi - expectedOffsetNoWriteLumi_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, global lumi";
      testsPass = false;
    }

    if (serviceTwo->nPreGlobalBeginLumi() != nGlobalBeginLumi ||
        serviceTwo->nPostGlobalBeginLumi() != nGlobalBeginLumi || serviceTwo->nPreGlobalEndLumi() != nGlobalBeginLumi ||
        serviceTwo->nPostGlobalEndLumi() != nGlobalBeginLumi ||
        serviceTwo->nPreModuleGlobalBeginLumi() != nGlobalBeginLumi * kNumberOfTestModules ||
        serviceTwo->nPostModuleGlobalBeginLumi() != nGlobalBeginLumi * kNumberOfTestModules ||
        serviceTwo->nPreModuleGlobalEndLumi() !=
            nGlobalBeginLumi * kNumberOfTestModules - expectedOffsetNoGlobalEndLumi_ ||
        serviceTwo->nPostModuleGlobalEndLumi() !=
            nGlobalBeginLumi * kNumberOfTestModules - expectedOffsetNoGlobalEndLumi_ ||
        serviceTwo->nPreGlobalWriteLumi() != nGlobalBeginLumi - expectedOffsetNoWriteLumi_ ||
        serviceTwo->nPostGlobalWriteLumi() != nGlobalBeginLumi - expectedOffsetNoWriteLumi_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceTwo, global lumi";
      testsPass = false;
    }

    if (serviceOne->nPreStreamBeginRun() != totalStreamBeginRun ||
        serviceOne->nPostStreamBeginRun() != totalStreamBeginRun ||
        serviceOne->nPreStreamEndRun() != totalStreamBeginRun ||
        serviceOne->nPostStreamEndRun() != totalStreamBeginRun ||
        serviceOne->nPreModuleStreamBeginRun() != totalStreamBeginRun * kNumberOfTestModules ||
        serviceOne->nPostModuleStreamBeginRun() != totalStreamBeginRun * kNumberOfTestModules ||
        serviceOne->nPreModuleStreamEndRun() !=
            totalStreamBeginRun * kNumberOfTestModules - expectedOffsetNoStreamEndRun_ ||
        serviceOne->nPostModuleStreamEndRun() !=
            totalStreamBeginRun * kNumberOfTestModules - expectedOffsetNoStreamEndRun_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, stream run";
      testsPass = false;
    }

    if (serviceTwo->nPreStreamBeginRun() != totalStreamBeginRun ||
        serviceTwo->nPostStreamBeginRun() != totalStreamBeginRun ||
        serviceTwo->nPreStreamEndRun() != totalStreamBeginRun ||
        serviceTwo->nPostStreamEndRun() != totalStreamBeginRun ||
        serviceTwo->nPreModuleStreamBeginRun() != totalStreamBeginRun * kNumberOfTestModules ||
        serviceTwo->nPostModuleStreamBeginRun() != totalStreamBeginRun * kNumberOfTestModules ||
        serviceTwo->nPreModuleStreamEndRun() !=
            totalStreamBeginRun * kNumberOfTestModules - expectedOffsetNoStreamEndRun_ ||
        serviceTwo->nPostModuleStreamEndRun() !=
            totalStreamBeginRun * kNumberOfTestModules - expectedOffsetNoStreamEndRun_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceTwo, stream run";
      testsPass = false;
    }

    unsigned int nGlobalBeginRun = nGlobalBeginRun_.load();

    if (serviceOne->nPreGlobalBeginRun() != nGlobalBeginRun || serviceOne->nPostGlobalBeginRun() != nGlobalBeginRun ||
        serviceOne->nPreGlobalEndRun() != nGlobalBeginRun || serviceOne->nPostGlobalEndRun() != nGlobalBeginRun ||
        serviceOne->nPreModuleGlobalBeginRun() != nGlobalBeginRun * kNumberOfTestModules ||
        serviceOne->nPostModuleGlobalBeginRun() != nGlobalBeginRun * kNumberOfTestModules ||
        serviceOne->nPreModuleGlobalEndRun() !=
            nGlobalBeginRun * kNumberOfTestModules - expectedOffsetNoGlobalEndRun_ ||
        serviceOne->nPostModuleGlobalEndRun() !=
            nGlobalBeginRun * kNumberOfTestModules - expectedOffsetNoGlobalEndRun_ ||
        serviceOne->nPreGlobalWriteRun() != nGlobalBeginRun - expectedOffsetNoWriteRun_ ||
        serviceOne->nPostGlobalWriteRun() != nGlobalBeginRun - expectedOffsetNoWriteRun_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceOne, global run";
      testsPass = false;
    }

    if (serviceTwo->nPreGlobalBeginRun() != nGlobalBeginRun || serviceTwo->nPostGlobalBeginRun() != nGlobalBeginRun ||
        serviceTwo->nPreGlobalEndRun() != nGlobalBeginRun || serviceTwo->nPostGlobalEndRun() != nGlobalBeginRun ||
        serviceTwo->nPreModuleGlobalBeginRun() != nGlobalBeginRun * kNumberOfTestModules ||
        serviceTwo->nPostModuleGlobalBeginRun() != nGlobalBeginRun * kNumberOfTestModules ||
        serviceTwo->nPreModuleGlobalEndRun() !=
            nGlobalBeginRun * kNumberOfTestModules - expectedOffsetNoGlobalEndRun_ ||
        serviceTwo->nPostModuleGlobalEndRun() !=
            nGlobalBeginRun * kNumberOfTestModules - expectedOffsetNoGlobalEndRun_ ||
        serviceTwo->nPreGlobalWriteRun() != nGlobalBeginRun - expectedOffsetNoWriteRun_ ||
        serviceTwo->nPostGlobalWriteRun() != nGlobalBeginRun - expectedOffsetNoWriteRun_) {
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "FAILED: Unexpected number of service transitions in TestServiceTwo, global run";
      testsPass = false;
    }

    if (verbose_) {
      edm::LogAbsolute("ExceptionThrowingProducer") << "nBeginJob_ = " << nBeginJob_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nBeginStream_ = " << nBeginStream_.load();
      edm::LogAbsolute("ExceptionThrowingProducer") << "nBeginProcessBlock_ = " << nBeginProcessBlock_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nEndProcessBlock_ = " << nEndProcessBlock_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nEndStream_ = " << nEndStream_.load();
      edm::LogAbsolute("ExceptionThrowingProducer") << "nEndJob_ = " << nEndJob_ << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer") << "nGlobalBeginLumi_ = " << nGlobalBeginLumi_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nGlobalEndLumi_ = " << nGlobalEndLumi_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nGlobalBeginRun_ = " << nGlobalBeginRun_;
      edm::LogAbsolute("ExceptionThrowingProducer") << "nGlobalEndRun_ = " << nGlobalEndRun_ << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPreBeginJob = " << serviceOne->nPreBeginJob();
      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPostBeginJob = " << serviceOne->nPostBeginJob();
      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPreEndJob = " << serviceOne->nPreEndJob();

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleBeginJob = " << serviceOne->nPreModuleBeginJob();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleBeginJob = " << serviceOne->nPostModuleBeginJob() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreBeginStream = " << serviceOne->nPreBeginStream();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostBeginStream = " << serviceOne->nPostBeginStream();
      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPreEndStream = " << serviceOne->nPreEndStream();
      edm::LogAbsolute("ExceptionThrowingProducer") << "serviceOne->nPostEndStream = " << serviceOne->nPostEndStream();

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleBeginStream = " << serviceOne->nPreModuleBeginStream();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleBeginStream = " << serviceOne->nPostModuleBeginStream();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleEndStream = " << serviceOne->nPreModuleEndStream();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleEndStream = " << serviceOne->nPostModuleEndStream() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreBeginProcessBlock = " << serviceOne->nPreBeginProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostBeginProcessBlock = " << serviceOne->nPostBeginProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreEndProcessBlock = " << serviceOne->nPreEndProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostEndProcessBlock = " << serviceOne->nPostEndProcessBlock();

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleBeginProcessBlock = " << serviceOne->nPreModuleBeginProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleBeginProcessBlock = " << serviceOne->nPostModuleBeginProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleEndProcessBlock = " << serviceOne->nPreModuleEndProcessBlock();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleEndProcessBlock = " << serviceOne->nPostModuleEndProcessBlock() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreStreamBeginLumi = " << serviceOne->nPreStreamBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostStreamBeginLumi = " << serviceOne->nPostStreamBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreStreamEndLumi = " << serviceOne->nPreStreamEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostStreamEndLumi = " << serviceOne->nPostStreamEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleStreamBeginLumi = " << serviceOne->nPreModuleStreamBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleStreamBeginLumi = " << serviceOne->nPostModuleStreamBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleStreamEndLumi = " << serviceOne->nPreModuleStreamEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleStreamEndLumi = " << serviceOne->nPostModuleStreamEndLumi() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalBeginLumi = " << serviceOne->nPreGlobalBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalBeginLumi = " << serviceOne->nPostGlobalBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalEndLumi = " << serviceOne->nPreGlobalEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalEndLumi = " << serviceOne->nPostGlobalEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleGlobalBeginLumi = " << serviceOne->nPreModuleGlobalBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleGlobalBeginLumi = " << serviceOne->nPostModuleGlobalBeginLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleGlobalEndLumi = " << serviceOne->nPreModuleGlobalEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleGlobalEndLumi = " << serviceOne->nPostModuleGlobalEndLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalWriteLumi = " << serviceOne->nPreGlobalWriteLumi();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalWriteLumi = " << serviceOne->nPostGlobalWriteLumi() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreStreamBeginRun = " << serviceOne->nPreStreamBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostStreamBeginRun = " << serviceOne->nPostStreamBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreStreamEndRun = " << serviceOne->nPreStreamEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostStreamEndRun = " << serviceOne->nPostStreamEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleStreamBeginRun = " << serviceOne->nPreModuleStreamBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleStreamBeginRun = " << serviceOne->nPostModuleStreamBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleStreamEndRun = " << serviceOne->nPreModuleStreamEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleStreamEndRun = " << serviceOne->nPostModuleStreamEndRun() << "\n";

      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalBeginRun = " << serviceOne->nPreGlobalBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalBeginRun = " << serviceOne->nPostGlobalBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalEndRun = " << serviceOne->nPreGlobalEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalEndRun = " << serviceOne->nPostGlobalEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleGlobalBeginRun = " << serviceOne->nPreModuleGlobalBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleGlobalBeginRun = " << serviceOne->nPostModuleGlobalBeginRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreModuleGlobalEndRun = " << serviceOne->nPreModuleGlobalEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostModuleGlobalEndRun = " << serviceOne->nPostModuleGlobalEndRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPreGlobalWriteRun = " << serviceOne->nPreGlobalWriteRun();
      edm::LogAbsolute("ExceptionThrowingProducer")
          << "serviceOne->nPostGlobalWriteRun = " << serviceOne->nPostGlobalWriteRun() << "\n";
    }

    if (testsPass) {
      edm::LogAbsolute("ExceptionThrowingProducer") << "All tests in ExceptionThrowingProducer PASSED";
    } else {
      edm::LogAbsolute("ExceptionThrowingProducer") << "At least one test in ExceptionThrowingProducer FAILED";
    }

    if (throwInEndJob_) {
      throw cms::Exception("IntentionalTestException")
          << "ExceptionThrowingProducer::endJob, module configured to throw during endJob";
    }
  }

  void ExceptionThrowingProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    edm::EventID invalidEventID;
    desc.addUntracked<bool>("verbose", false);
    desc.addUntracked<edm::EventID>("eventIDThrowOnEvent", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnGlobalBeginRun", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnGlobalBeginLumi", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnGlobalEndRun", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnGlobalEndLumi", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnStreamBeginRun", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnStreamBeginLumi", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnStreamEndRun", invalidEventID);
    desc.addUntracked<edm::EventID>("eventIDThrowOnStreamEndLumi", invalidEventID);

    desc.addUntracked<bool>("throwInBeginJob", false);
    desc.addUntracked<bool>("throwInBeginStream", false);
    desc.addUntracked<bool>("throwInBeginProcessBlock", false);
    desc.addUntracked<bool>("throwInEndProcessBlock", false);
    desc.addUntracked<bool>("throwInEndStream", false);
    desc.addUntracked<bool>("throwInEndJob", false);

    desc.addUntracked<unsigned int>("expectedNBeginJob", 1);
    desc.addUntracked<unsigned int>("expectedNBeginStream", 4);
    desc.addUntracked<unsigned int>("expectedNBeginProcessBlock", 1);
    desc.addUntracked<unsigned int>("expectedNEndProcessBlock", 1);
    desc.addUntracked<unsigned int>("expectedNEndStream", 4);
    desc.addUntracked<unsigned int>("expectedNEndJob", 1);
    desc.addUntracked<bool>("expectNoRunsProcessed", false);
    desc.addUntracked<unsigned int>("expectedOffsetNoEndJob", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoEndStream", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoEndProcessBlock", 0);
    desc.addUntracked<unsigned int>("expectedStreamBeginLumi", kUnset);
    desc.addUntracked<unsigned int>("expectedOffsetNoStreamEndLumi", 0);
    desc.addUntracked<unsigned int>("expectedGlobalBeginLumi", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoGlobalEndLumi", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoWriteLumi", 0);

    desc.addUntracked<unsigned int>("expectedStreamBeginRun", kUnset);
    desc.addUntracked<unsigned int>("expectedOffsetNoStreamEndRun", 0);
    desc.addUntracked<unsigned int>("expectedGlobalBeginRun", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoGlobalEndRun", 0);
    desc.addUntracked<unsigned int>("expectedOffsetNoWriteRun", 0);

    descriptions.addDefault(desc);
  }

}  // namespace edmtest
using edmtest::ExceptionThrowingProducer;
DEFINE_FWK_MODULE(ExceptionThrowingProducer);
