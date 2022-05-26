// Generated from C:/Users/Filippo/Projects/TertiaryStructuresComparator/src/main/resources\TertiaryStructureBonds.g4 by ANTLR 4.10.1

	package it.unicam.cs.bdslab.tersaling;

import org.antlr.v4.runtime.tree.ParseTreeVisitor;

/**
 * This interface defines a complete generic visitor for a parse tree produced
 * by {@link TertiaryStructureBondsParser}.
 *
 * @param <T> The return type of the visit operation. Use {@link Void} for
 * operations with no return type.
 */
public interface TertiaryStructureBondsVisitor<T> extends ParseTreeVisitor<T> {
	/**
	 * Visit a parse tree produced by the {@code bondsContinue}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitBondsContinue(TertiaryStructureBondsParser.BondsContinueContext ctx);
	/**
	 * Visit a parse tree produced by the {@code bondsEnd}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitBondsEnd(TertiaryStructureBondsParser.BondsEndContext ctx);
	/**
	 * Visit a parse tree produced by {@link TertiaryStructureBondsParser#bond}.
	 * @param ctx the parse tree
	 * @return the visitor result
	 */
	T visitBond(TertiaryStructureBondsParser.BondContext ctx);
}