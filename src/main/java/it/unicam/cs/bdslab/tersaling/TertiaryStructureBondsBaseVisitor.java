// Generated from C:/Users/Filippo/Projects/TertiaryStructuresComparator/src/main/resources\TertiaryStructureBonds.g4 by ANTLR 4.10.1

	package it.unicam.cs.bdslab.tersaling;

import org.antlr.v4.runtime.tree.AbstractParseTreeVisitor;

/**
 * This class provides an empty implementation of {@link TertiaryStructureBondsVisitor},
 * which can be extended to create a visitor which only needs to handle a subset
 * of the available methods.
 *
 * @param <T> The return type of the visit operation. Use {@link Void} for
 * operations with no return type.
 */
public class TertiaryStructureBondsBaseVisitor<T> extends AbstractParseTreeVisitor<T> implements TertiaryStructureBondsVisitor<T> {
	/**
	 * {@inheritDoc}
	 *
	 * <p>The default implementation returns the result of calling
	 * {@link #visitChildren} on {@code ctx}.</p>
	 */
	@Override public T visitBondsContinue(TertiaryStructureBondsParser.BondsContinueContext ctx) { return visitChildren(ctx); }
	/**
	 * {@inheritDoc}
	 *
	 * <p>The default implementation returns the result of calling
	 * {@link #visitChildren} on {@code ctx}.</p>
	 */
	@Override public T visitBondsEnd(TertiaryStructureBondsParser.BondsEndContext ctx) { return visitChildren(ctx); }
	/**
	 * {@inheritDoc}
	 *
	 * <p>The default implementation returns the result of calling
	 * {@link #visitChildren} on {@code ctx}.</p>
	 */
	@Override public T visitBond(TertiaryStructureBondsParser.BondContext ctx) { return visitChildren(ctx); }
}