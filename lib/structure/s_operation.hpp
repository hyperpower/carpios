#ifndef _S_OPERATION_HPP_
#define _S_OPERATION_HPP_

#include <structure/s_data.hpp>
#include "s_define.hpp"
#include "s_grid.hpp"
#include "s_vector.hpp"
#include "s_data.hpp"

namespace structure {

template<St DIM>
class Operation_ {
public:
	static const St Dim = DIM;
	typedef Equation_<DIM> Base;
	typedef Equation_<DIM> Equation;
	typedef Grid_<DIM> Grid;
	typedef typename Grid::Poi Poi;
	typedef std::shared_ptr<Grid> spGrid;
	typedef Index_<Dim> Index;
	typedef Ijk_<Dim> Ijk;

	typedef Scalar_<DIM> Scalar;
	typedef std::shared_ptr<Scalar> spScalar;

	typedef VectorCenter_<DIM> VecC;
	typedef std::shared_ptr<VecC> spVecC;

	typedef VectorFace_<DIM> VecF;
	typedef std::shared_ptr<VecF> spVecF;

	typedef std::function<Vt(Vt, Vt, Vt, Vt)> Function;

public:
	static void Set(spScalar spdata, Function fun) {
		Scalar& cdata = *(spdata);
		spGrid grid = cdata.get_grid();
		for (typename Grid::Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Vt val = fun(0, grid->C_(_X_, IJK.i()), grid->C_(_Y_, IJK.j()),
					grid->C_(_Z_, IJK.k()));
			cdata.VAL(IJK) = val;
		}
	}

	static spVecF Grad(spScalar spcdata) {	// Return grad(f) (Result: vFace)
		const Scalar& cdata = *(spcdata);
		spGrid grid = cdata.get_grid();
		// build spVecf
		spScalar asp[] = { nullptr, nullptr, nullptr };
		for (St d = 0; d < DIM; ++d) {
			asp[d] = spScalar(new Scalar(grid));
		}
		spVecF spvface = spVecF(new VecF(asp[0], asp[1], asp[2]));
		VecF& vface = *spvface;
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Vt vc = cdata(idx);
				Index idxp = idx.p(d);
				vface[d](idx) = (cdata(idxp) - vc)
						/ (grid->c_(d, idxp) - grid->c_(d, idx));
				if (idx[d] == 0) {
					Index idxm = idx.m(d);
					vface[d](idxm) = (cdata(idxm) - vc)
							/ (grid->c_(d, idxm) - grid->c_(d, idx));
				}
			}
		}
		return spvface;
	}
	static spVecC GradCenter(spScalar spcdata) {// Return grad(f) (Result: vCell)
		const Scalar& cdata = *(spcdata);
		spGrid grid = cdata.get_grid();
		// build spVecf
		spScalar asp[] = { nullptr, nullptr, nullptr };
		for (St d = 0; d < DIM; ++d) {
			asp[d] = spScalar(new Scalar(grid));
		}
		spVecC spvcenter = spVecC(new VecC(asp[0], asp[1], asp[2]));
		VecC& vcenter = *spvcenter;
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Vt vc = cdata(idx);
				Index idxp = idx.p(d);
				Index idxm = idx.m(d);
				Vt gradfp = (cdata(idxp) - vc)
						/ (grid->c_(d, idxp) - grid->c_(d, idx));
				Vt gradfm = (vc - cdata(idxm))
						/ (grid->c_(d, idx) - grid->c_(d, idxm));
				vcenter[d](idx) = (gradfp + gradfm) * 0.5;
			}
		}
		return spvcenter;
	}

	static spScalar Div(spVecF spvf) {
		spGrid grid = spvf->get_grid();
		const VecF& vface = *spvf;
		spScalar spscalar = spScalar(new Scalar(grid));
		Scalar& scalar = *spscalar;
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Vt div = 0;
			Index idx = ijk.current();
			//std::cout<< idx <<std::endl;
			for (St d = 0; d < DIM; ++d) {
				Index idxm = idx.m(d);
				Vt vpf = vface[d](idx);
				Vt vmf = vface[d](idxm);
				Vt s = grid->s_(d, idx[d]);
				//std::cout<< "d" <<d<< " vpf "<<vpf << " vmf "<<vmf<<std::endl;
				div += (vpf - vmf) / s;
			}
			scalar(idx) = div;
		}
		return spscalar;
	}
	static spScalar Div(spScalar sps) {
		spGrid grid = sps->get_grid();
		const Scalar& scalar = *sps;
		spScalar spscalar = spScalar(new Scalar(grid));
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Vt div = 0;
			Index idx = ijk.current();
			for (St d = 0; d < DIM; ++d) {
				Index idxp = idx.p(d);
				Index idxm = idx.m(d);
				Vt vp = scalar(idxp);
				Vt vm = scalar(idxm);
				Vt dis = grid->c_(d, idxp[d]) - grid->c_(d, idxm[d]);
				div += (vp - vm) / dis;
			}
			(*spscalar)(idx) = div;
		}
		return spscalar;
	}

	static void InterpolateF2C(VecF& vface, VecC& vcenter) {// Return grad(f) (Result: vCell)
		spGrid grid = vcenter[0].get_grid();
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Index idxm = idx.m(d);
				vcenter[d](idx) = (vface[d](idx) + vface[d](idxm)) * 0.5;
			}
		}
	}
	static spVecC InterpolateF2C(spVecF vface) {// Return grad(f) (Result: vCell)
		spGrid grid = vface->get_grid();
		// build spVecf
		spScalar asp[] = { nullptr, nullptr, nullptr };
		for (St d = 0; d < DIM; ++d) {
			asp[d] = spScalar(new Scalar(grid));
		}
		spVecC spvcenter = spVecC(new VecC(asp[0], asp[1], asp[2]));
		VecC& vcenter = *spvcenter;
		// ------------
		InterpolateF2C(*vface, vcenter);
		return spvcenter;
	}

	static void InterpolateC2F(spVecC spcdata, spVecF spvface) {// Return grad(f) (Result: vCell)
		VecC& vcenter = *spcdata;
		VecF& vface = *spvface;
		spGrid grid = vcenter[0].get_grid();
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Index idxp = idx.p(d);
				Scalar& cdata = vcenter[d];
				Vt vc = cdata(idx);
				Vt hs = cdata.hs_(d, idx);
				Vt hsp = cdata.hs_(d, idxp);
				vface[d](idx) = (cdata(idxp) * hs + vc * hsp) / (hs + hsp);
				if (idx[d] == 0) {
					Index idxm = idx.m(d);
					Vt hsm = cdata.hs_(d, idxm);
					vface[d](idxm) = (cdata(idxm) * hs + vc * hsm) / (hs + hsm);
				}
			}
		}
	}

protected:

public:
	static Vt InterpolateCoordinate(spScalar spcdata, Vt x, Vt y = 0.0, Vt z =
			0.0) {
		const Scalar& cdata = *(spcdata);
		spGrid grid = cdata.get_grid();
		// 1 the point should be in the range of grid
		Poi p(x, y, z);
		if (grid->is_in_on(p)) {
			// 1d -> 2
			// 2d -> 4
			// 3d -> 8
			Index aindex[2];
			for (St d = 0; d < DIM; ++d) {
				aindex[0][d] = grid->find_close_idx_m(d, p[d]);
				aindex[1][d] = grid->find_close_idx_p(d, p[d]);
			}

			Poi pmm = grid->c(aindex[0]);
			Poi ppp = grid->c(aindex[1]);

			Arr aav(std::pow(2, DIM));
			if (DIM == 1) {
				aav[0] = cdata.val(aindex[0]);
				aav[1] = cdata.val(aindex[1]);
			}
			if (DIM == 2) {
				aav[0] = cdata.val(aindex[0]);
				aav[1] = cdata.val(aindex[0][0], aindex[1][1]);
				aav[2] = cdata.val(aindex[1]);
				aav[3] = cdata.val(aindex[1][0], aindex[0][1]);
				//std::cout << aav[0] << " ";
				//std::cout << aav[1] << " ";
				//std::cout << aav[2] << " ";
				//std::cout << aav[3] << "\n";
			}
			if (DIM == 3) {
				aav[0] = cdata.val(aindex[0]);
				aav[1] = cdata.val(aindex[0][0], aindex[1][1], aindex[0][2]);
				aav[2] = cdata.val(aindex[1][0], aindex[1][1], aindex[0][2]);
				aav[3] = cdata.val(aindex[1][0], aindex[0][1], aindex[0][2]);
				aav[4] = cdata.val(aindex[0][0], aindex[1][1], aindex[1][2]);
				aav[5] = cdata.val(aindex[0][0], aindex[1][1], aindex[1][2]);
				aav[6] = cdata.val(aindex[1]);
				aav[7] = cdata.val(aindex[1][0], aindex[0][1], aindex[1][2]);
			}
			//std::cout << aindex[0] << " ";
			//std::cout << aindex[1] << "\n";
			//std::cout << pmm << " ";
			//std::cout << ppp << "\n";
			return carpio::Interp_<Vt, Vt>::Linear_all(DIM, p, pmm, ppp, aav);

		} else {
			return 0;
		}
	}

	static spVecF InterpolateC2F(spScalar spcdata) {// Return grad(f) (Result: vCell)
		const Scalar& cdata = *(spcdata);
		spGrid grid = cdata.get_grid();
		// build spVecf
		spScalar asp[] = { nullptr, nullptr, nullptr };
		for (St d = 0; d < DIM; ++d) {
			asp[d] = spScalar(new Scalar(grid));
		}
		spVecF spvface = spVecF(new VecF(asp[0], asp[1], asp[2]));
		VecF& vface = *spvface;
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Vt vc = cdata(idx);
				Index idxp = idx.p(d);
				Vt hs = grid->hs_(d, idx[d]);
				Vt hsp = grid->hs_(d, idxp[d]);
				vface[d](idx) = (cdata(idxp) * hs + vc * hsp) / (hs + hsp);
				if (idx[d] == 0) {
					Index idxm = idx.m(d);
					Vt hsm = grid->hs_(d, idxm[d]);
					vface[d](idxm) = (cdata(idxm) * hs + vc * hsm) / (hs + hsm);
				}
			}
		}
		return spvface;
	}

	static void InterpolateC2F(spScalar spcdata, spVecF spvface) {// Return grad(f) (Result: vCell)
		const Scalar& cdata = *(spcdata);
		spGrid grid = cdata.get_grid();
		// build spVecf
		//spScalar asp[] = { nullptr, nullptr, nullptr };
		//for (St d = 0; d < DIM; ++d) {
		//	asp[d] = spScalar(new Scalar(grid));
		//}
		//spVecF spvface = spVecF(new VecF(asp[0], asp[1], asp[2]));
		VecF& vface = *spvface;
		// ------------
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			for (St d = 0; d < DIM; ++d) {
				Index idx = ijk.current();
				Vt vc = cdata(idx);
				Index idxp = idx.p(d);
				Vt hs = grid->hs_(d, idx[d]);
				Vt hsp = grid->hs_(d, idxp[d]);
				vface[d](idx) = (cdata(idxp) * hs + vc * hsp) / (hs + hsp);
				if (idx[d] == 0) {
					Index idxm = idx.m(d);
					Vt hsm = grid->hs_(d, idxm[d]);
					vface[d](idxm) = (cdata(idxm) * hs + vc * hsm) / (hs + hsm);
				}
			}
		}
		//return spvface;
	}

	static spScalar VdotNabla(spVecF spvf, spScalar sps,
			const std::string& name) {
		if (name == "UPWIND1") {
			return VdotNabla_upwind1(spvf, sps);
		}
		SHOULD_NOT_REACH;
		return nullptr;
	}

	static spScalar VdotNabla_upwind1(const VecF& vface, const Scalar& cdata) {
		spGrid grid = cdata.get_grid();
		spScalar spscalar = spScalar(new Scalar(grid));
		// build spVecf
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Index idx = ijk.current();
			Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
			for (St d = 0; d < DIM; ++d) {
				Index idxp = idx.p(d);
				Index idxm = idx.m(d);
				Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
				if (veo > 0.0) {
					v_dot_dfd[d] = veo * (cdata(idx) - cdata(idxm))
							/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
				} else if (veo < 0.0) {
					v_dot_dfd[d] = veo * (cdata(idxp) - cdata(idx))
							/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
				}
			}
			// ** summary ** //
			Vt v_dot_nabla_f = 0.0;
			for (St d = 0; d < DIM; ++d) {
				v_dot_nabla_f += v_dot_dfd[d];
			}
			(*spscalar)(idx) = v_dot_nabla_f;
		}
		return spscalar;
	}

	static spScalar VdotNabla_upwind1(spVecF spvf, spScalar sps) {
		Scalar& cdata = *(sps);
		spGrid grid = cdata.get_grid();
		spScalar spscalar = spScalar(new Scalar(grid));
		// build spVecf
		VecF& vface = *spvf;
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Index idx = ijk.current();
			Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
			for (St d = 0; d < DIM; ++d) {
				Index idxp = idx.p(d);
				Index idxm = idx.m(d);
				Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
				if (veo > 0.0) {
					v_dot_dfd[d] = veo * (cdata(idx) - cdata(idxm))
							/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
				} else if (veo < 0.0) {
					v_dot_dfd[d] = veo * (cdata(idxp) - cdata(idx))
							/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
				}
			}
			// ** summary ** //
			Vt v_dot_nabla_f = 0.0;
			for (St d = 0; d < DIM; ++d) {
				v_dot_nabla_f += v_dot_dfd[d];
			}
			(*spscalar)(idx) = v_dot_nabla_f;
		}
		return spscalar;
	}

	static Vt VdotNabla_upwind1(spVecF spvf, spScalar sps, Index idx) {
		Scalar& cdata = *(sps);
		spScalar spscalar = spScalar(new Scalar(cdata.get_grid()));
		// build spVecf
		VecF& vface = *spvf;
		//
		Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);
			Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
			if (veo > 0.0) {
				v_dot_dfd[d] = veo * (cdata(idx) - cdata(idxm))
						/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
			} else if (veo < 0.0) {
				v_dot_dfd[d] = veo * (cdata(idxp) - cdata(idx))
						/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
			}
		}
		// ** summary ** //
		Vt v_dot_nabla_f = 0.0;
		for (St d = 0; d < DIM; ++d) {
			v_dot_nabla_f += v_dot_dfd[d];
		}
		return v_dot_nabla_f;
	}

	static Vt VdotNabla_upwind1(const VecF& vface, const Scalar& cdata,
			Index idx) {
		//
		Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);
			Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
			if (veo > 0.0) {
				v_dot_dfd[d] = veo * (cdata(idx) - cdata(idxm))
						/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
			} else if (veo < 0.0) {
				v_dot_dfd[d] = veo * (cdata(idxp) - cdata(idx))
						/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
			}
		}
		// ** summary ** //
		Vt v_dot_nabla_f = 0.0;
		for (St d = 0; d < DIM; ++d) {
			v_dot_nabla_f += v_dot_dfd[d];
		}
		return v_dot_nabla_f;
	}

	static Vt VdotNabla_upwind1(const VecF& vface, const Scalar& cdata,
			Index idx, St d) {
		//
		Vt v_dot_dfd = 0.0;
		ASSERT(d < DIM);
		Index idxp = idx.p(d);
		Index idxm = idx.m(d);
		Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
		if (veo > 0.0) {
			v_dot_dfd = veo * (cdata(idx) - cdata(idxm))
					/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
		} else if (veo < 0.0) {
			v_dot_dfd = veo * (cdata(idxp) - cdata(idx))
					/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
		}
		return v_dot_dfd;
	}

	static Vt VdotNabla_center(const VecF& vface, const Scalar& cdata,
			Index idx) {
		//
		Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);
			Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
			v_dot_dfd[d] = veo * (cdata(idxp) - cdata(idxm))
					/ (cdata.c_(d, idxp) - cdata.c_(d, idxm));
		}
		// ** summary ** //
		Vt v_dot_nabla_f = 0.0;
		for (St d = 0; d < DIM; ++d) {
			v_dot_nabla_f += v_dot_dfd[d];
		}
		return v_dot_nabla_f;
	}

	static Vt VdotNabla_center4(const VecF& vface, const Scalar& cdata,
			Index idx) {
		Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxpp = idxp.p(d);
			Index idxm = idx.m(d);
			Index idxmm = idxm.m(d);
			Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
			Vt dfdx_inn = (cdata(idxp) - cdata(idxm))
					/ (cdata.c_(d, idxp) - cdata.c_(d, idxm));
			Vt dfdx_out = (cdata(idxpp) - cdata(idxmm))
					/ (cdata.c_(d, idxpp) - cdata.c_(d, idxmm));
			v_dot_dfd[d] = (4.0 * dfdx_inn - dfdx_out) / 3.0;
		}
		// ** summary ** //
		Vt v_dot_nabla_f = 0.0;
		for (St d = 0; d < DIM; ++d) {
			v_dot_nabla_f += v_dot_dfd[d];
		}
		return v_dot_nabla_f;
	}

	static Vt _rCD(const Scalar& cdata, St d, const Index& U, const Index& C,
			const Index& D) {
		Vt sU = cdata.s_(d, U);
		Vt sC = cdata.s_(d, C);
		Vt sD = cdata.s_(d, D);
		Vt vU = cdata(U);
		Vt vC = cdata(C);
		Vt vD = cdata(D);
		return (vC - vU) * (sD + sC) / (vD - vC + SMALL) / (sC + sU);
	}

	static Vt _RCD(const Scalar& cdata, St d, const Index& C, const Index& D) {
		Vt sC = cdata.s_(d, C);
		Vt sD = cdata.s_(d, D);
		return (sD + sC) / sC;
	}

	static Vt limiter_VanLeer(Vt r, Vt R) {
		return 0.5 * R * (r + std::abs(r)) / (R - 1 + r);
	}

	static Vt limiter_superbee(Vt r, Vt R) {
		return std::max(std::max(0.0, std::min(R * r, 1.0)), std::min(r, R));
	}

	static Vt limiter_WAHYD(Vt r, Vt R) {
		if (r <= 1.0) {
			return 0.5 * R * (r + std::abs(r)) / (R - 1 + r);
		} else {
			return std::min((r + R * r * std::abs(r)) / (R + r * r), R);
		}
	}

	//  A review on TVD schemes and a refined flux-limiter
	//  for steady-state calculations
	//  Di Zhang, Chunbo Jiang, Dongfang Liang, Liang Cheng
	//  Journal of Computational Physics 302 (2015) 114–154
	//
	// k-scheme
	//            1 + k           1 - k
	// limiter = ------- r(CD) + -------
	//              2               2
	//          d(phi)/dx (CU)
	// r(CD) = ----------------
	//          d(phi)/dx (DC)
	//
	// SOU                      k = -1  (Second order upwind   upwind2)
	// Fromm                    k = 0
	// CUI                      k = 1/3
	// QUICK                    k = 1/2
	// CDS                      k = 1   (Center Difference Scheme,  center)

	static Vt limiter_QUICK(Vt r, Vt R) {
		return 0.75 * r + 0.25;
	}
	static Vt limiter_CUI(Vt r, Vt R) {
		return 2.0 / 3.0 * r + 1.0 / 3.0;
	}

	typedef Vt (*Limiter)(Vt, Vt);

	// Improved total variation diminishing schemes for advection
	// simulation on arbitrary grids
	// J. Hou * ,† , F. Simons and R. Hinkelmann
	// Int. J. Numer. Meth. Fluids 2012; 70:359–382
	static Vt VdotNabla_TVD(const VecF& vface, const Scalar& cdata, Index idx,
			Limiter lim) {
		//
		Vt v_dot_dfd[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);
			Vt veo = (vface[d](idx) + vface[d](idxm)) * 0.5;
			Index Cm,Um,Dm,Cp,Up,Dp;
			Vt fp,fm,R,r;
			if (veo >= 0.0) {
				// fm ------------------
				Dm = idx;
				Cm = idxm;
				Um = idxm.m(d);
				// fp ------------------
				Dp = idxp;
				Cp = idx;
				Up = idxm;
			} else if (veo < 0.0) {
				// fm ------------------
				Dm = idxm;
				Cm = idx;
				Um = idxp;
				// fp ------------------
				Dp = idx;
				Cp = idxp;
				Up = idxp.p(d);
			}
			// p -------------------------
			R = _RCD(cdata, d, Cp, Dp);
			r = _rCD(cdata, d, Up, Cp, Dp);
			fp = cdata(Cp) + lim(r, R) / R * (cdata(Dp) - cdata(Cp));
			// m ---------------------
			R = _RCD(cdata, d, Cm, Dm);
			r = _rCD(cdata, d, Um, Cm, Dm);
			fm = cdata(Cm) + lim(r, R) / R * (cdata(Dm) - cdata(Cm));
			//
			v_dot_dfd[d] = (vface[d](idx) * fp - vface[d](idxm) * fm)
					/ cdata.s_(d, idx);
		}
		// ** summary ** //
		Vt v_dot_nabla_f = 0.0;
		for (St d = 0; d < DIM; ++d) {
			v_dot_nabla_f += v_dot_dfd[d];
		}
		return v_dot_nabla_f;
	}

	static Vt limiter_upwind2(Vt r, Vt R) {
		return r;
	}

	static Vt VdotNabla_upwind2(const VecF& vface, const Scalar& cdata,
			Index idx) {
		// there are no limiter
		return VdotNabla_TVD(vface, cdata, idx, limiter_upwind2);
	}

	static Vt VdotNabla_QUICK(const VecF& vface, const Scalar& cdata,
			Index idx) {
		return VdotNabla_TVD(vface, cdata, idx, limiter_QUICK);
	}

	static Vt VdotNabla_CUI(const VecF& vface, const Scalar& cdata, Index idx) {
		return VdotNabla_TVD(vface, cdata, idx, limiter_CUI);
	}

	static Vt VdotNabla_TVD_VanLeer(const VecF& vface, const Scalar& cdata,
			Index idx) {
		return VdotNabla_TVD(vface, cdata, idx, limiter_VanLeer);
	}
	static Vt VdotNabla_TVD_superbee(const VecF& vface, const Scalar& cdata,
			Index idx) {
		return VdotNabla_TVD(vface, cdata, idx, limiter_superbee);
	}
	static Vt VdotNabla_TVD_WAHYD(const VecF& vface, const Scalar& cdata,
			Index idx) {
		return VdotNabla_TVD(vface, cdata, idx, limiter_WAHYD);
	}

	static Vt NablaMuNabla(const Scalar& cmu, const Scalar& cdata, Index idx) {
		//Scalar& cdata = *(sps);
		//Scalar& cmu = *(spmu);
		// build spVecf
		//
		Vt tdata[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);

			Vt dfdx_p,dfdx_m;
			Vt mufp,mufm;
			mufp = (cdata.hs_(d, idx) * cmu(idxp)
					+ cdata.hs_(d, idxp) * cmu(idx))
					/ (cdata.hs_(d, idxp) + cdata.hs_(d, idx));
			mufm = (cdata.hs_(d, idx) * cmu(idxm)
					+ cdata.hs_(d, idxm) * cmu(idx))
					/ (cdata.hs_(d, idxm) + cdata.hs_(d, idx));

			dfdx_p = mufp * (cdata(idxp) - cdata(idx))
					/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
			dfdx_m = mufm * (cdata(idx) - cdata(idxm))
					/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
			Vt d2fdx2 = (dfdx_p - dfdx_m) / cdata.s_(d, idx);
			tdata[d] = d2fdx2; //tmp.SetData(m,d2fdx2);
		}
		// ** summary ** //
		Vt res = 0.0;
		for (St d = 0; d < DIM; ++d) {
			res += tdata[d];
		}
		return res;
	}

	static Vt NablaMuNabla(Vt mu, const Scalar& cdata, Index idx) {
		// Scalar& cdata = *(sps);
		// build spVecf
		//
		Vt tdata[] = { 0.0, 0.0, 0.0 };
		for (St d = 0; d < DIM; ++d) {
			Index idxp = idx.p(d);
			Index idxm = idx.m(d);

			Vt dfdx_p,dfdx_m;
			dfdx_m = mu * (cdata(idx) - cdata(idxm))
					/ (cdata.c_(d, idx) - cdata.c_(d, idxm));
			dfdx_p = mu * (cdata(idxp) - cdata(idx))
					/ (cdata.c_(d, idxp) - cdata.c_(d, idx));
			Vt d2fdx2 = (dfdx_p - dfdx_m) / cdata.s_(d, idx);

			tdata[d] = d2fdx2; //tmp.SetData(m,d2fdx2);
		}
		// ** summary ** //
		Vt res = 0.0;
		for (St d = 0; d < DIM; ++d) {
			res += tdata[d];
		}
		return res;
	}

	static void NablaMuNabla(Vt mu, const Scalar& scalar, Scalar& res) {
		Grid& grid = *(scalar.get_grid());
		for (Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index idx = ijk.current();
			Vt val = NablaMuNabla(mu, scalar, idx);
			res(idx) = val;
		}
	}
	static spScalar NablaMuNabla(Vt mu, const Scalar& scalar) {
		spGrid grid = scalar.get_grid();
		spScalar spres= spScalar(new Scalar(grid));
		NablaMuNabla(mu, scalar, *spres);
		return spres;
	}

	static void NablaMuNabla(const Scalar& mu, const Scalar& scalar,
			Scalar& res) {
		Grid& grid = *(scalar.get_grid());
		for (Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index idx = ijk.current();
			Vt val = NablaMuNabla(mu, scalar, idx);
			res(idx) = val;
		}
	}

	static spScalar NablaMuNabla(const Scalar& mu, const Scalar& scalar) {
			spGrid grid = scalar.get_grid();
			spScalar spres= spScalar(new Scalar(grid));
			NablaMuNabla(mu, scalar, *spres);
			return spres;
	}

	static void Divide(spScalar sps, Vt value) {
		ASSERT(value != 0.0);
		spGrid grid = sps->get_grid();
		Scalar& scalar = *sps;
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& idx = ijk.current();
			scalar(idx) = scalar(idx) / value;
		}
	}

	static void Multiply(spScalar sps, Vt value) {
		Scalar& scalar = *sps;
		Multiply(scalar, value);
	}

	static void Multiply(Scalar& scalar, Vt value) {
		spGrid grid = scalar.get_grid();
		for (Ijk ijk = grid->begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& idx = ijk.current();
			scalar(idx) = scalar(idx) * value;
		}
	}

	static void Multiply(spVecF sps, Vt value) {
		spGrid grid = sps->get_grid();
		VecF& vfa = *sps;
		// ------------
		for (St d = 0; d < DIM; ++d) {
			Multiply(vfa[d], value);
		}
	}

	static void Minus(spScalar spa, spScalar spb) {
		spGrid grid = spa->get_grid();
		Scalar& scalara = *spa;
		Scalar& scalarb = *spb;
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			scalara.VAL(IDX) = scalara.VAL(IDX) - scalarb.VAL(IDX);
		}
	}
	static void Minus(Scalar& a, const Scalar& b) {
		spGrid grid = a.get_grid();
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			a.VAL(IDX) = a.VAL(IDX) - b.VAL(IDX);
		}
	}

	static void Add(Scalar& a, const Scalar& b) {
		spGrid grid = a.get_grid();
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			a.VAL(IDX) = a.VAL(IDX) + b.VAL(IDX);
		}
	}

	static void Add(Scalar& a, const Vt& b) {
		spGrid grid = a.get_grid();
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			a.VAL(IDX) = a.VAL(IDX) + b;
		}
	}

	static spScalar new_Minus(spScalar spa, spScalar spb) {
		spGrid grid = spa->get_grid();
		spScalar res = spScalar(new Scalar(grid));
		Scalar& scalara = *spa;
		Scalar& scalarb = *spb;
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			res->VAL(IDX) = scalara.VAL(IDX) - scalarb.VAL(IDX);
		}
		return res;
	}
	// spa = spa - spb
	static void Minus(spVecF spa, spVecF spb) {
		spGrid grid = spa->get_grid();
		VecF& vfa = *spa;
		VecF& vfb = *spb;
		// ------------
		for (St d = 0; d < DIM; ++d) {
			Minus(vfa[d], vfb[d]);
		}
	}
	static void Minus(spVecC spa, spVecC spb) {
		spGrid grid = spa->get_grid();
		VecC& vca = *spa;
		VecC& vcb = *spb;
		// ------------
		for (St d = 0; d < DIM; ++d) {
			Minus(vca[d], vcb[d]);
		}
	}
	// b = a
	static void Copy(Scalar& a, Scalar& b) {
		spGrid grid = a.get_grid();
		for (Ijk IJK = grid->begin_IJK(); !IJK.is_end(); ++IJK) {
			Index& IDX = IJK.current();
			b.VAL(IDX) = a.VAL(IDX);
		}
	}
	static void Copy(spScalar spa, spScalar& spb) {
		Scalar& a = *spa;
		Scalar& b = *spb;
		Copy(a, b);
	}
	static void Copy(spVecF spa, spVecF spb) {
		VecF& vfa = *spa;
		VecF& vfb = *spb;
		// ------------
		for (St d = 0; d < DIM; ++d) {
			Copy(vfa[d], vfb[d]);
		}
	}
	static void Copy(spVecC spa, spVecC spb) {
		VecC& va = *spa;
		VecC& vb = *spb;
		// ------------
		for (St d = 0; d < DIM; ++d) {
			Copy(va[d], vb[d]);
		}
	}

	static Vt TimeRestrict_CourantNumber(Vt courantNumberAllowed, VecC& vc) {
		Vt timeStepRestricted;
		Grid& grid = *(vc.get_grid());
		Vt maxMag = vc.max_magnitude();
		Vt minSize = grid.min_size();
		timeStepRestricted = courantNumberAllowed * minSize / maxMag;
		return timeStepRestricted;
	}

	static Vt TimeRestrict_DiffusionNumber(Vt diffusionNumberAllowed,
			const Grid& grid, Vt diffusionCoefficient) {
		Vt timeStepRestricted;
		Vt denominator = 0.0;
		for (St d = 0; d < Dim; d++) {
			Vt dmin = grid.min_size(d);
			denominator += 1.0 / (dmin * dmin);
		}
		timeStepRestricted = diffusionNumberAllowed / diffusionCoefficient
				/ denominator;
		return timeStepRestricted;
	}

	static Vt Error1(const Scalar& exa, const Scalar& val) {
		Grid& grid = (*exa.get_grid());
		Vt norm = 0;
		Vt svol = 0;
		for (Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& idx = ijk.current();
			Vt vol = grid.volume(idx);
			Vt err = val(idx) - exa(idx);
			norm += (std::abs(err) * vol);
			svol += vol;
		}
		return norm / svol;
	}

	static Vt Error2(const Scalar& exa, const Scalar& val) {
		Grid& grid = (*exa.get_grid());
		Vt norm = 0;
		Vt svol = 0;
		for (Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& idx = ijk.current();
			Vt vol = grid.volume(idx);
			Vt err = val(idx) - exa(idx);
			norm += (err * err * vol);
			svol += vol;
		}
		return std::sqrt(norm) / svol;
	}

	static Vt Errori(const Scalar& exa, const Scalar& val) {
		Grid& grid = (*exa.get_grid());
		Vt normi = 0;
		for (Ijk ijk = grid.begin_ijk(); !ijk.is_end(); ++ijk) {
			Index& idx = ijk.current();
			Vt err = std::abs(val(idx) - exa(idx));
			if (err > normi) {
				normi = err;
			}
		}
		return normi;
	}

};

}

#endif
