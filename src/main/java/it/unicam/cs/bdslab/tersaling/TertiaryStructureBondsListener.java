// Generated from C:/Users/Filippo/Projects/TertiaryStructuresComparator/src/main/resources\TertiaryStructureBonds.g4 by ANTLR 4.10.1

	package it.unicam.cs.bdslab.tersaling;

import org.antlr.v4.runtime.tree.ParseTreeListener;

/**
 * This interface defines a complete listener for a parse tree produced by
 * {@link TertiaryStructureBondsParser}.
 */
public interface TertiaryStructureBondsListener extends ParseTreeListener {
	/**
	 * Enter a parse tree produced by the {@code bondsContinue}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 */
	void enterBondsContinue(TertiaryStructureBondsParser.BondsContinueContext ctx);
	/**
	 * Exit a parse tree produced by the {@code bondsContinue}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 */
	void exitBondsContinue(TertiaryStructureBondsParser.BondsContinueContext ctx);
	/**
	 * Enter a parse tree produced by the {@code bondsEnd}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 */
	void enterBondsEnd(TertiaryStructureBondsParser.BondsEndContext ctx);
	/**
	 * Exit a parse tree produced by the {@code bondsEnd}
	 * labeled alternative in {@link TertiaryStructureBondsParser#bonds}.
	 * @param ctx the parse tree
	 */
	void exitBondsEnd(TertiaryStructureBondsParser.BondsEndContext ctx);
	/**
	 * Enter a parse tree produced by {@link TertiaryStructureBondsParser#bond}.
	 * @param ctx the parse tree
	 */
	void enterBond(TertiaryStructureBondsParser.BondContext ctx);
	/**
	 * Exit a parse tree produced by {@link TertiaryStructureBondsParser#bond}.
	 * @param ctx the parse tree
	 */
	void exitBond(TertiaryStructureBondsParser.BondContext ctx);
}