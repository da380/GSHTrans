// The parallel loops: every index exactly once, the sequential exception
// whatever the team and schedule, and a reduction with the same bits on any
// team.

#include <gtest/gtest.h>

#include <GSHTrans/Core.hpp>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

using namespace GSHTrans;
using Int = std::ptrdiff_t;

namespace {

const auto teams = std::vector<Execution>{
    Execution::Sequential(), Execution::Parallel(2), Execution::Parallel(3),
    Execution::Parallel(5),  Execution::Parallel(8), Execution::Parallel(16)};

std::string Describe(Execution policy) {
  return policy.IsParallel() ? std::to_string(policy.Threads()) + " threads"
                             : "sequential";
}

// Throws at every index in `failures`, with the index as the message.
struct Failing {
  std::vector<Int> failures;
  void operator()(Int i) const {
    for (auto f : failures) {
      if (i == f) throw std::runtime_error(std::to_string(i));
    }
  }
};

std::string MessageOf(auto&& call) {
  try {
    call();
  } catch (const std::runtime_error& error) {
    return error.what();
  }
  return "no exception";
}

}  // namespace

TEST(Parallel, EveryIndexExactlyOnce) {
  for (auto n : {Int{0}, Int{1}, Int{7}, Int{1000}}) {
    for (auto policy : teams) {
      auto seen = std::vector<std::atomic<int>>(static_cast<std::size_t>(n));
      ForEach(policy, n, [&](Int i) { seen[static_cast<std::size_t>(i)]++; });
      for (auto& s : seen)
        EXPECT_EQ(s.load(), 1) << n << ", " << Describe(policy);

      auto tasks = std::vector<std::atomic<int>>(static_cast<std::size_t>(n));
      ForEachTask(policy, n,
                  [&](Int i) { tasks[static_cast<std::size_t>(i)]++; });
      for (auto& s : tasks)
        EXPECT_EQ(s.load(), 1) << n << ", " << Describe(policy);

      auto blocks = std::vector<std::atomic<int>>(static_cast<std::size_t>(n));
      ForEachBlock(policy, n, [&](Int begin, Int end) {
        for (auto i = begin; i < end; ++i)
          blocks[static_cast<std::size_t>(i)]++;
      });
      for (auto& s : blocks)
        EXPECT_EQ(s.load(), 1) << n << ", " << Describe(policy);
    }
  }
}

// The exception a sequential loop would stop at -- the lowest failing index --
// whichever thread fails first in time, for every form.
TEST(Parallel, RethrowsTheLowestFailingIndex) {
  constexpr auto n = Int{1000};
  const auto body = Failing{{900, 37, 500}};
  for (auto policy : teams) {
    EXPECT_EQ(MessageOf([&] { ForEach(policy, n, body); }), "37")
        << Describe(policy);
    EXPECT_EQ(MessageOf([&] { ForEachTask(policy, n, body); }), "37")
        << Describe(policy);
    EXPECT_EQ(MessageOf([&] {
                ForEachBlock(policy, n, [&](Int begin, Int end) {
                  for (auto i = begin; i < end; ++i) body(i);
                });
              }),
              "37")
        << Describe(policy);
    EXPECT_EQ(MessageOf([&] {
                Reduce(
                    policy, n, 0.0,
                    [&](Int i) {
                      body(i);
                      return 1.0;
                    },
                    [](double a, double b) { return a + b; });
              }),
              "37")
        << Describe(policy);
  }
}

// What makes the lowest failure findable: nothing below it is skipped.
TEST(Parallel, RunsEveryIndexBelowAFailure) {
  constexpr auto n = Int{2000};
  constexpr auto failure = Int{1500};
  for (auto policy : teams) {
    auto seen = std::vector<std::atomic<int>>(static_cast<std::size_t>(n));
    const auto body = [&](Int i) {
      seen[static_cast<std::size_t>(i)]++;
      if (i == failure) throw std::runtime_error("failed");
    };
    EXPECT_THROW(ForEachTask(policy, n, body), std::runtime_error);
    for (auto i = Int{0}; i < failure; ++i) {
      ASSERT_EQ(seen[static_cast<std::size_t>(i)].load(), 1)
          << "index " << i << ", " << Describe(policy);
    }
  }
}

// Terms spread over thirty orders of magnitude, so that any change of
// grouping changes the bits.
TEST(Parallel, ReduceHasTheSameBitsOnAnyTeam) {
  const auto term = [](Int i) {
    return std::sin(0.37 * static_cast<double>(i)) *
           std::pow(10.0, static_cast<double>(i % 31) - 15.0);
  };
  for (auto n : {Int{1}, Int{100}, Int{255}, Int{256}, Int{100003}}) {
    const auto reference = Sum(Execution::Sequential(), n, 0.0, term);
    for (auto policy : teams) {
      EXPECT_EQ(Sum(policy, n, 0.0, term), reference)
          << n << " terms, " << Describe(policy);
    }
  }
  EXPECT_EQ(Sum(Execution::Parallel(4), Int{0}, 2.5, term), 2.5);
}

TEST(Parallel, ReduceCombinesInitOnceAndFirst) {
  // A non-commutative combine shows the order: concatenation.
  const auto digit = [](Int i) { return std::to_string(i % 10); };
  const auto concat = [](const std::string& a, const std::string& b) {
    return a + b;
  };
  auto expected = std::string{"<"};
  for (auto i = Int{0}; i < 600; ++i) expected += digit(i);
  for (auto policy : teams) {
    EXPECT_EQ(Reduce(policy, Int{600}, std::string{"<"}, digit, concat),
              expected)
        << Describe(policy);
  }
}

// One level threads: inside an active region a loop runs on the thread that
// called it.
TEST(Parallel, RunsSequentiallyInsideARegion) {
  auto ok = std::atomic<bool>{true};
#pragma omp parallel num_threads(4)
  {
    const auto caller = OpenMP::ThreadNumber();
    ForEach(Execution::Parallel(4), Int{64}, [&](Int) {
      if (OpenMP::ThreadNumber() != caller) ok = false;
    });
  }
  EXPECT_TRUE(ok.load());
}
