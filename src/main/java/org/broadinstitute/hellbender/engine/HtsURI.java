package org.broadinstitute.hellbender.engine;

import java.net.URI;
import java.nio.file.Path;

/**
 * Interface for htsjdk input/output URIs.
 *
 * // TODO: should this have a getInputStream (only if getPath() == true)
 * // TODO: should there be a method that returns a "file" path
 */
public interface HtsURI {

    /**
     * Return true if this {code HtsURI} can be resolved to an NIO Path. If true, {@code #toPath()} can be safely called.
     *
     * There are cases where a valid URI with a valid scheme backed by an installed NIO File System
     * still can't be turned into a Path, i.e., the following specifies an invalid authority "namenode":
     *
     *  hdfs://namenode/to/file
     *
     * The current implementation returns false for these cases (toPath will fail, getInvalidPathReason
     * returns the reason code).
     */
    boolean isPath();

    /**
     * @return true if this URI has a scheme that has an installed NIO file system provider. This does not guarantee
     * the URI can be converted into a Path, since the URI can be syntactically valid, and specify a valid file
     * system provider, but still fail to be semantically meaningful.
     */
    boolean isNIO();

    /**
     * Resolve this HtsURI to an NIO Path. Can be safely called only if {@code #isPath()} returns true.
     */
    Path toPath();

    /**
     * Return a string message describing why this URI cannot be converted to a Path ({@code #isPath()} returns false).
     * @return Message explaining toPath failure reason, since it can fail for various reasons.
     */
    String getToPathFailureReason();

    /**
     * Get a {@code java.net.URI} object for this {@code HtsURI}. Will not be null.
     * @return The {@code URI} object for this HtsURI.
     */
    URI getURI();

    /**
     * Returns the string from which this {@code HtsURI} was originally created. This string may differ from the normalized
     * string returned from a Path that has been object resolved from this HtsURI.
     * //TODO: should this return a file:// scheme for URIs with no other scheme
     * @return string from which this URI as originally created. Will not be null.
     */
    String getURIString();

    /**
     * Return the scheme for this HtsURI. For file URIs (URIs that have no explicit scheme), this will return
     * the scheme "file".
     * @return the scheme String or this URI, if any. May be null.
     */
    default String getScheme() {
        return getURI().getScheme();
    }
}
