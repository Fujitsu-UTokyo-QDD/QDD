
#include "common.h"
#include "dd.h"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <random>

#define M_PI 3.14159265358979323846

static vEdge apply_gate(QubitCount count, vEdge v, GateMatrix g, Qubit target) {
  auto gate = makeGate(count, g, target);
  return mv_multiply(gate, v);
}
static vEdge apply_gate_with_control(QubitCount count, vEdge v, GateMatrix g,
                                     Qubit target, Qubit control) {
  Controls controls;
  controls.emplace(Control{control, Control::Type::pos});
  auto gate = makeGate(count, g, target, controls);
  return mv_multiply(gate, v);
}

static std::uint64_t gcd(std::uint64_t a, std::uint64_t b) {
  while (a != 0) {
    const std::uint64_t c = a;
    a = b % a;
    b = c;
  }
  return b;
}

static std::uint64_t modpow(std::uint64_t base, std::uint64_t exp,
                            std::uint64_t modulus) {
  base %= modulus;
  std::uint64_t result = 1ULL;
  while (exp > 0) {
    if ((exp & 1ULL) > 0) {
      result = (result * base) % modulus;
    }
    base = (base * base) % modulus;
    exp >>= 1ULL;
  }
  return result;
}

static double cosine(double fac, double div) {
  return std::cos((M_PI * fac) / div);
}

static double sine(double fac, double div) {
  return std::sin((M_PI * fac) / div);
}

std::pair<std::uint32_t, std::uint32_t>
postProcessing(const std::string &sample, std::size_t requiredBits,
               std::size_t compositeN, std::size_t coprimeA) {
  std::ostream log{nullptr};

  std::uint64_t res = 0;

  log << "measurement: ";
  for (std::uint32_t i = 0; i < 2 * requiredBits; i++) {
    log << sample.at(2 * requiredBits - 1 - i);
    res = (res << 1U) + (sample.at(requiredBits + i) == '1' ? 1 : 0);
  }

  log << " = " << res << "\n";

  if (res == 0) {
    log << "Factorization failed (measured 0)!\n";
    return {0, 0};
  }
  std::vector<std::uint64_t> cf;
  auto denom = 1ULL << (2 * requiredBits);
  const auto oldDenom = denom;
  const auto oldRes = res;

  log << "Continued fraction expansion of " << res << "/" << denom << " = "
      << std::flush;

  while (res != 0) {
    cf.push_back(denom / res);
    const auto tmp = denom % res;
    denom = res;
    res = tmp;
  }

  for (const auto i : cf) {
    log << i << " ";
  }
  log << "\n";

  for (std::uint32_t i = 0; i < cf.size(); i++) {
    // determine candidate
    std::uint64_t denominator = cf[i];
    std::uint64_t numerator = 1;

    for (std::int32_t j = static_cast<std::int32_t>(i) - 1; j >= 0; j--) {
      const auto tmp =
          numerator + cf[static_cast<std::size_t>(j)] * denominator;
      numerator = denominator;
      denominator = tmp;
    }

    log << "  Candidate " << numerator << "/" << denominator << ": ";

    if (denominator > compositeN) {
      log << " denominator too large (greater than " << compositeN
          << ")!\nFactorization failed!\n";
      return {0, 0};
    }

    const double delta =
        static_cast<double>(oldRes) / static_cast<double>(oldDenom) -
        static_cast<double>(numerator) / static_cast<double>(denominator);
    if (std::abs(delta) >= 1.0 / (2.0 * static_cast<double>(oldDenom))) {
      log << "delta is too big (" << delta << ")\n";
      continue;
    }

    std::uint64_t fact = 1;
    while (denominator * fact < compositeN &&
           modpow(coprimeA, denominator * fact, compositeN) != 1) {
      fact++;
    }

    if (modpow(coprimeA, denominator * fact, compositeN) != 1) {
      log << "failed\n";
      continue;
    }

    log << "found period: " << denominator << " * " << fact << " = "
        << (denominator * fact) << "\n";

    if (((denominator * fact) & 1U) > 0) {
      log << "Factorization failed (period is odd)!\n";
      return {0, 0};
    }

    auto f1 = modpow(coprimeA, (denominator * fact) / 2, compositeN);
    auto f2 = (f1 + 1) % compositeN;
    f1 = (f1 == 0) ? compositeN - 1 : f1 - 1;
    f1 = gcd(f1, compositeN);
    f2 = gcd(f2, compositeN);

    if (f1 == 1ULL || f2 == 1ULL) {
      log << "Factorization failed: found trivial factors " << f1 << " and "
          << f2 << "\n";
      return {0, 0};
    }
    log << "Factorization succeeded! Non-trivial factors are: \n"
        << "  -- gcd(" << compositeN << "^(" << (denominator * fact) << "/2)-1"
        << "," << compositeN << ") = " << f1 << "\n"
        << "  -- gcd(" << compositeN << "^(" << (denominator * fact) << "/2)+1"
        << "," << compositeN << ") = " << f2 << "\n";
    return {f1, f2};
  }
  return {0, 0};
}

mEdge limitTo(std::uint64_t a, std::size_t requiredBits) {
  std::array<mEdge, 4> edges{mEdge::zero, mEdge::zero, mEdge::zero,
                             mEdge::zero};

  if ((a & 1U) > 0) {
    edges[0] = edges[3] = mEdge::one;
  } else {
    edges[0] = mEdge::one;
  }
  mEdge f = makeMEdge(0, edges);

  edges[0] = edges[1] = edges[2] = edges[3] = mEdge::zero;

  for (std::uint32_t p = 1; p < requiredBits + 1; p++) {
    if (((a >> p) & 1U) > 0) {
      edges[0] = makeIdent(p - 1);
      edges[3] = f;
    } else {
      edges[0] = f;
    }
    f = makeMEdge(static_cast<Qubit>(p), edges);
    edges[3] = mEdge::zero;
  }

  return f;
}

mEdge addConst(std::uint64_t a, std::size_t requiredBits) {
  mEdge f = mEdge::one;
  std::array<mEdge, 4> edges{mEdge::zero, mEdge::zero, mEdge::zero,
                             mEdge::zero};

  std::uint32_t p = 0;
  while (((a >> p) & 1U) == 0U) {
    edges[0] = f;
    edges[3] = f;
    f = makeMEdge(static_cast<Qubit>(p), edges);
    p++;
  }

  mEdge right;
  mEdge left;

  edges[0] = edges[1] = edges[2] = edges[3] = mEdge::zero;
  edges[2] = f;
  left = makeMEdge(static_cast<Qubit>(p), edges);
  edges[2] = mEdge::zero;
  edges[1] = f;
  right = makeMEdge(static_cast<Qubit>(p), edges);
  p++;

  mEdge newLeft;
  mEdge newRight;
  for (; p < requiredBits; p++) {
    edges[0] = edges[1] = edges[2] = edges[3] = mEdge::zero;
    if (((a >> p) & 1U) != 0U) {
      edges[2] = left;
      newLeft = makeMEdge(static_cast<Qubit>(p), edges);
      edges[2] = mEdge::zero;
      edges[0] = right;
      edges[1] = left;
      edges[3] = right;
      newRight = makeMEdge(static_cast<Qubit>(p), edges);
    } else {
      edges[1] = right;
      newRight = makeMEdge(static_cast<Qubit>(p), edges);
      edges[1] = mEdge::zero;
      edges[0] = left;
      edges[2] = right;
      edges[3] = left;
      newLeft = makeMEdge(static_cast<Qubit>(p), edges);
    }
    left = newLeft;
    right = newRight;
  }

  edges[0] = left;
  edges[1] = right;
  edges[2] = right;
  edges[3] = left;

  return makeMEdge(static_cast<Qubit>(p), edges);
}

mEdge addConstMod(std::uint64_t a, std::size_t compositeN,
                  std::size_t requiredBits) {
  const auto f = addConst(a, requiredBits);
  const auto f2 = addConst(compositeN, requiredBits);
  const auto f3 = limitTo(compositeN - 1, requiredBits);

  auto f4 = limitTo(compositeN - 1 - a, requiredBits);
  f4.w = f4.w.neg();

  const auto diff2 = mm_add(f3, f4);

  f4.w = f4.w.neg();

  const auto simEdgeMultiply = mm_multiply(conjugateTranspose(f2), f);
  const auto simEdgeResult =
      mm_add(mm_multiply(f, f4), mm_multiply(simEdgeMultiply, diff2));

  return simEdgeResult.n->children[0];
}

vEdge uAEmulate(vEdge rootEdge, std::size_t nQubits, std::size_t requiredBits,
                std::uint64_t a, std::int32_t q, std::size_t compositeN) {
  const mEdge limit = makeIdent(static_cast<Qubit>(requiredBits));

  mEdge f = mEdge::one;
  std::array<mEdge, 4> edges{mEdge::zero, mEdge::zero, mEdge::zero,
                             mEdge::zero};

  for (std::uint32_t p = 0; p < requiredBits; ++p) {
    edges[0] = f;
    edges[1] = f;
    f = makeMEdge(static_cast<Qubit>(p), edges);
  }

  // TODO: limitTo?

  f = mm_multiply(f, limit);

  edges[1] = mEdge::zero;

  std::uint64_t t = a;

  for (std::uint32_t i = 0; i < requiredBits; ++i) {
    mEdge active = mEdge::one;
    for (std::uint32_t p = 0; p < requiredBits; ++p) {
      if (p == i) {
        edges[3] = active;
        edges[0] = mEdge::zero;
      } else {
        edges[0] = edges[3] = active;
      }
      active = makeMEdge(static_cast<Qubit>(p), edges);
    }

    active.w = Complex(-1.0, 0.0);
    const mEdge passive = mm_multiply(f, mm_add(limit, active));
    active.w = Complex(1.0, 0.0);
    active = mm_multiply(f, active);

    const mEdge tmp = addConstMod(t, compositeN, requiredBits);
    active = mm_multiply(tmp, active);

    f = mm_add(active, passive);

    t = (2 * t) % compositeN;
  }

  mEdge e = f;

  for (auto i = static_cast<std::int32_t>(2 * requiredBits - 1); i >= 0; --i) {
    if (i == q) {
      edges[1] = edges[2] = mEdge::zero;
      edges[0] = mEdge::one;
      edges[3] = e;
      e = makeMEdge(
          static_cast<Qubit>(nQubits - 1 - static_cast<std::size_t>(i)), edges);
    } else {
      edges[1] = edges[2] = mEdge::zero;
      edges[0] = edges[3] = e;
      e = makeMEdge(
          static_cast<Qubit>(nQubits - 1 - static_cast<std::size_t>(i)), edges);
    }
  }

  const vEdge tmp = mv_multiply(e, rootEdge);

  return tmp;
}

std::pair<Complex, std::string> getPathOfLeastResistance(std::size_t nQubits,
                                                         vEdge rootEdge) {
  double epsilon = std::numeric_limits<double>::epsilon();
  if (std::abs(rootEdge.w.mag2() - 1.0) > epsilon) {
    if (rootEdge.w.isApproximatelyZero()) {
      throw std::runtime_error(
          "Numerical instabilities led to a 0-vector! Abort simulation!");
    }
    std::cerr << "WARNING in PoLR: numerical instability occurred during "
                 "simulation: |alpha|^2 + |beta|^2 - 1 = "
              << 1.0 - rootEdge.w.mag2() << ", but should be 1!\n";
  }

  std::string result(nQubits, '0');
  auto pathValue = rootEdge.w;
  vEdge cur = rootEdge;
  for (Qubit i = rootEdge.n->v + 1; i-- > 0;) {
    pathValue = pathValue * cur.w;
    double p0 = cur.n->children.at(0).w.mag2();
    double p1 = cur.n->children.at(1).w.mag2();
    double tmp = p0 + p1;

    if (std::abs(tmp - 1.0) > epsilon) {
      throw std::runtime_error("Added probabilities differ from 1 by " +
                               std::to_string(std::abs(tmp - 1.0)));
    }
    p0 /= tmp; // TODO: should p1 be normalized as well?

    if (p0 >= p1) {
      cur = cur.n->children.at(0);
    } else {
      result[static_cast<std::size_t>(cur.n->v)] = '1';
      cur = cur.n->children.at(1);
    }
  }

  return {pathValue, std::string{result.rbegin(), result.rend()}};
}
std::map<std::string, std::size_t>
shor_simulate(std::size_t compositeN, std::size_t coprimeA, std::size_t shots) {
  std::mt19937_64 mt;
  auto requiredBits =
      static_cast<std::size_t>(std::ceil(std::log2(compositeN)));
  QubitCount nQubits = static_cast<QubitCount>(3 * requiredBits);
  auto rootEdge = makeZeroState(nQubits);
  // Initialize qubits
  // TODO: other init method where the initial value can be chosen
  auto applyGate = [nQubits](vEdge v, GateMatrix g, Qubit target) {
    return apply_gate(nQubits, v, g, target);
  };
  auto applyGateWithControl = [nQubits](vEdge v, GateMatrix g, Qubit target,
                                        Qubit control) {
    return apply_gate_with_control(nQubits, v, g, target, control);
  };
  rootEdge = applyGate(rootEdge, Xmat, 0);

  if (coprimeA != 0 && gcd(coprimeA, compositeN) != 1) {

    coprimeA = 0;
  }

  if (coprimeA == 0) {
    std::uniform_int_distribution<std::size_t> distribution(
        1, compositeN - 1); // range is inclusive
    do {
      coprimeA = distribution(mt);
    } while (gcd(coprimeA, compositeN) != 1 || coprimeA == 1);
  }

  std::vector<std::uint64_t> as;
  as.resize(2 * requiredBits);
  assert(as.size() == 2 * requiredBits); // it's quite easy to get the vector
                                         // initialization wrong
  as[2 * requiredBits - 1] = coprimeA;
  std::uint64_t newA = coprimeA;
  for (auto i = static_cast<std::int64_t>(2 * requiredBits - 2); i >= 0; i--) {
    newA = newA * newA;
    newA = newA % compositeN;
    as[static_cast<std::size_t>(i)] = newA;
  }

  for (std::uint32_t i = 0; i < 2 * requiredBits; i++) {
    rootEdge = applyGate(rootEdge, Hmat, static_cast<Qubit>((nQubits - 1) - i));
  }
  const auto mod = static_cast<std::int32_t>(
      std::ceil(2.0 * static_cast<double>(requiredBits) /
                6.0)); // log_0.9(0.5) is about 6
  const auto t1 = std::chrono::steady_clock::now();

  for (std::uint32_t i = 0; i < 2 * requiredBits; i++) {
    rootEdge = uAEmulate(rootEdge, nQubits, requiredBits, as[i],
                         static_cast<std::int32_t>(i), compositeN);
  }

  // EXACT QFT
  for (std::int32_t i = 0; i < static_cast<std::int32_t>(2 * requiredBits);
       i++) {

    double q = 2;

    for (std::int32_t j = i - 1; j >= 0; j--) {
      double qR = cosine(1, -q);
      double qI = sine(1, -q);
      const GateMatrix qm{1, 0, 0, {qR, qI}};
      rootEdge = applyGateWithControl(
          rootEdge, qm,
          static_cast<Qubit>(nQubits - 1 - static_cast<std::size_t>(i)),
          static_cast<Qubit>(nQubits - 1 - static_cast<std::size_t>(j)));
      q *= 2;
    }

    rootEdge = applyGate(
        rootEdge, Hmat,
        static_cast<Qubit>(nQubits - 1 - static_cast<std::size_t>(i)));
  }

  {
    std::string sampleReversed = measureAll(rootEdge, false, mt);
    const std::string sample{sampleReversed.rbegin(), sampleReversed.rend()};
    auto simFactors =
        postProcessing(sample, requiredBits, compositeN, coprimeA);
    std::string simResult;
    if (simFactors.first != 0 && simFactors.second != 0) {
      simResult = std::string("SUCCESS(") + std::to_string(simFactors.first) +
                  "*" + std::to_string(simFactors.second) + ")";
    } else {
      simResult = "FAILURE";
    }
    std::cout << simResult << std::endl;
  }

  // the path of least resistance result (does not involve randomness)
  {
    const std::pair<Complex, std::string> polrPair =
        getPathOfLeastResistance(static_cast<std::size_t>(nQubits), rootEdge);
    std::clog << polrPair.first << " " << polrPair.second << "\n";
    std::string polrReversed = polrPair.second;
    const std::string polr = {polrReversed.rbegin(), polrReversed.rend()};
    auto polrFactors = postProcessing(polr, requiredBits, compositeN, coprimeA);
    std::string polrResult;
    if (polrFactors.first != 0 && polrFactors.second != 0) {
      polrResult = std::string("SUCCESS(") + std::to_string(polrFactors.first) +
                   "*" + std::to_string(polrFactors.second) + ")";
    } else {
      polrResult = "FAILURE";
    }
    std::cout << polrResult << std::endl;
  }

  return {};
}