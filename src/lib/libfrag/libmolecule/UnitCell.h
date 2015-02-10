/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#ifndef SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_
#define SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_

#include "Molecule.h"
namespace psi{
namespace LibMolecule{

class UnitCell: public Molecule{
   private:
      ///The angles of the unit cell, in Radians
      double angles_[3];
      ///The lengths of the sides of the unit cell, in a.u.
      double sides_[3];
      ///The fractional coordinate to cart transformation matrix
      double Frac2Cart_[9];
      ///The cart w fractional transform
      double Cart2Frac_[9];
   protected:
      void SetTrans();
   public:
      void SetAngles(const double alpha,const double beta,
                     const double gamma,const bool IsDegree=true);
      void SetSides(const double a, const double b, const double c,
                    const bool IsBohr=false);
      UnitCell(boost::shared_ptr<Molecule> Mol,const bool IsFrac=true);
};

}}



#endif /* SRC_LIB_LIBFRAG_LIBMOLECULE_UNITCELL_H_ */